import sys
from pathlib import Path
from unittest.mock import MagicMock
import numpy as np
import pytest
sys.path.append('../')
from pyRichards import richards
from pyRichards import mssubsurface


exchange_smc_regional_units = mssubsurface.exchange_smc_regional_units
exchange_dz_regional_units = mssubsurface.exchange_dz_regional_units

class FakeBus:
    def __init__(self):
        self.messages = {}
        self.sent_log = []

    def preload(self, source, dest, tag, payload):
        print(
            f"[FakeBus.preload] source={source} dest={dest} tag={tag} shape={np.array(payload).shape}",
            flush=True,
        )
        self.messages.setdefault((source, dest, tag), []).append(np.array(payload, copy=True))

    def push(self, source, dest, tag, payload):
        payload_copy = np.array(payload, copy=True)
        print(
            f"[FakeBus.push] source={source} dest={dest} tag={tag} shape={payload_copy.shape}",
            flush=True,
        )
        self.sent_log.append((source, dest, tag, payload_copy))
        self.messages.setdefault((source, dest, tag), []).append(payload_copy)

    def pop(self, source, dest, tag):
        print(f"[FakeBus.pop] source={source} dest={dest} tag={tag}", flush=True)
        key = (source, dest, tag)
        if key not in self.messages or not self.messages[key]:
            raise AssertionError(
                f"No queued message for source={source}, dest={dest}, tag={tag}"
            )
        payload = self.messages[key].pop(0)
        print(
            f"[FakeBus.pop] delivered shape={payload.shape} remaining={len(self.messages[key])}",
            flush=True,
        )
        return payload


class FakeRequest:
    def __init__(self, kind, bus, rank, peer, tag, buffer_ref=None):
        self.kind = kind
        self.bus = bus
        self.rank = rank
        self.peer = peer
        self.tag = tag
        self.buffer_ref = buffer_ref

    def complete(self):
        print(
            f"[FakeRequest.complete] kind={self.kind} rank={self.rank} peer={self.peer} tag={self.tag}",
            flush=True,
        )
        if self.kind == "recv":
            incoming = self.bus.pop(source=self.peer, dest=self.rank, tag=self.tag)
            self.buffer_ref[:, :] = incoming
            print(
                f"[FakeRequest.complete] recv buffer filled shape={self.buffer_ref.shape}",
                flush=True,
            )


class FakeComm:
    def __init__(self, bus, rank, allgather_remote=None):
        self.bus = bus
        self.rank = rank
        self.allgather_remote = list(allgather_remote) if allgather_remote is not None else []

    def Isend(self, payload, dest, tag):
        print(f"[FakeComm.Isend] rank={self.rank} -> dest={dest} tag={tag}", flush=True)
        self.bus.push(source=self.rank, dest=dest, tag=tag, payload=payload)
        return FakeRequest(kind="send", bus=self.bus, rank=self.rank, peer=dest, tag=tag)

    def Irecv(self, buffer_ref, source, tag):
        print(
            f"[FakeComm.Irecv] rank={self.rank} <- source={source} tag={tag} expected_shape={buffer_ref.shape}",
            flush=True,
        )
        return FakeRequest(
            kind="recv",
            bus=self.bus,
            rank=self.rank,
            peer=source,
            tag=tag,
            buffer_ref=buffer_ref,
        )

    def Barrier(self):
        print(f"[FakeComm.Barrier] rank={self.rank}", flush=True)
        return

    def allgather(self, value):
        print(f"[FakeComm.allgather] rank={self.rank} local_len={len(value)}", flush=True)
        gathered = [value]
        if self.allgather_remote:
            gathered.append(self.allgather_remote.pop(0))
        return gathered


class DummySubsurface:
    def __init__(self, comm, th_gw, risfu_mapping, reg_ids, nsoil, dz_gw=None):
        self.comm = comm
        self.th_gw = th_gw
        self.dz_gw = dz_gw if dz_gw is not None else np.array(th_gw, copy=True)
        self.risfu_mapping = risfu_mapping
        self.reg_ids = reg_ids
        self.nsoil = nsoil

        nrows = sum(len(reg_ids[cid]) for cid in reg_ids)
        self.reg_theta_gw = np.zeros((nrows, nsoil), dtype=float)
        self.reg_dz_gw = np.zeros((nrows, nsoil), dtype=float)


class DummyHB:
    def __init__(self, subsurface, cid_rank_mapping):
        self.mssubsurface = subsurface
        self.cid_rank_mapping = cid_rank_mapping


def test_exchange_smc_regional_units_external_exchange_and_assembly(monkeypatch):
    print("\n[test_exchange_smc_regional_units_external_exchange_and_assembly] start", flush=True)
    bus = FakeBus()

    # Emulate rank 0 processing cid 1 with remote dependency on cid 2
    rank = 0
    cids = [1]
    cid_rank_mapping = {1: 0, 2: 1}

    # cid 1 has two local units with 2 soil layers
    th_gw_local = np.array(
        [
            [10.0, 100.0],
            [11.0, 110.0],
        ]
    )

    # Global mapping from prior risfu_connections_regional step
    # - cid1 needs rows [0] from cid2
    # - cid2 needs rows [1] from cid1
    risfu_mapping = {
        1: {
            1: np.array([0, 1], dtype=int),
            2: np.array([0], dtype=int),
        },
        2: {
            1: np.array([1], dtype=int),
            2: np.array([0, 1], dtype=int),
        },
    }

    # Assembly order: local cid 1 first (2 rows), then remote cid 2 (1 row)
    reg_ids = {
        1: np.array([1, 2], dtype=int),
        2: np.array([1], dtype=int),
    }

    smc_tag_offset = 200000
    remote_send_plan = [(1, 0, smc_tag_offset + 2001, 1, 2, str(th_gw_local.dtype))]
    remote_recv_plan = [(0, 1, smc_tag_offset + 1002, 1, 2, str(th_gw_local.dtype))]
    comm = FakeComm(bus=bus, rank=rank, allgather_remote=[remote_send_plan, remote_recv_plan])
    subsurface = DummySubsurface(
        comm=comm,
        th_gw=th_gw_local,
        risfu_mapping=risfu_mapping,
        reg_ids=reg_ids,
        nsoil=2,
    )

    HBdb = {
        1: DummyHB(subsurface=subsurface, cid_rank_mapping=cid_rank_mapping),
    }

    # Incoming payload from rank 1 (cid 2 -> cid 1) uses SMC tag namespace offset.
    remote_payload = np.array([[20.0, 200.0]])
    print(f"[test] preloading remote payload: {remote_payload}", flush=True)
    bus.preload(source=1, dest=0, tag=smc_tag_offset + 2001, payload=remote_payload)

    def fake_waitall(requests):
        print(f"[fake_waitall] n_requests={len(requests)}", flush=True)
        for req in requests:
            req.complete()

    class _FakeMPI:
        class Request:
            @staticmethod
            def Waitall(requests):
                fake_waitall(requests)

    monkeypatch.setattr(mssubsurface, "MPI", _FakeMPI)

    exchange_smc_regional_units(cids=cids, rank=rank, HBdb=HBdb)
    print(f"[test] sent_log entries={len(bus.sent_log)}", flush=True)

    # Verify send behavior: cid 2 requested local row [1], sent with SMC-offset tag 201002
    sends = [
        entry
        for entry in bus.sent_log
        if entry[0] == 0 and entry[1] == 1 and entry[2] == smc_tag_offset + 1002
    ]
    print(f"[test] matching sends count={len(sends)}", flush=True)
    assert len(sends) == 1
    print(f"[test] send payload found={sends[0][3]}", flush=True)
    np.testing.assert_allclose(sends[0][3], np.array([[11.0, 110.0]]))

    # Verify final assembled regional theta matrix
    expected_reg_theta = np.array(
        [
            [10.0, 100.0],
            [11.0, 110.0],
            [20.0, 200.0],
        ]
    )
    print(f"[test] expected reg_theta_gw=\n{expected_reg_theta}", flush=True)
    print(f"[test] actual reg_theta_gw=\n{HBdb[1].mssubsurface.reg_theta_gw}", flush=True)
    np.testing.assert_allclose(HBdb[1].mssubsurface.reg_theta_gw, expected_reg_theta)
    print("[test_exchange_smc_regional_units_external_exchange_and_assembly] done", flush=True)


def test_exchange_dz_regional_units_external_exchange_and_assembly(monkeypatch):
    print("\n[test_exchange_dz_regional_units_external_exchange_and_assembly] start", flush=True)
    bus = FakeBus()

    rank = 0
    cids = [1]
    cid_rank_mapping = {1: 0, 2: 1}

    dz_gw_local = np.array(
        [
            [0.5, 1.5],
            [0.6, 1.6],
        ]
    )

    risfu_mapping = {
        1: {
            1: np.array([0, 1], dtype=int),
            2: np.array([0], dtype=int),
        },
        2: {
            1: np.array([1], dtype=int),
            2: np.array([0, 1], dtype=int),
        },
    }

    reg_ids = {
        1: np.array([1, 2], dtype=int),
        2: np.array([1], dtype=int),
    }

    dz_tag_offset = 300000
    remote_send_plan = [(1, 0, dz_tag_offset + 2001, 1, 2, str(dz_gw_local.dtype))]
    remote_recv_plan = [(0, 1, dz_tag_offset + 1002, 1, 2, str(dz_gw_local.dtype))]
    comm = FakeComm(bus=bus, rank=rank, allgather_remote=[remote_send_plan, remote_recv_plan])
    subsurface = DummySubsurface(
        comm=comm,
        th_gw=np.array(dz_gw_local, copy=True),
        dz_gw=dz_gw_local,
        risfu_mapping=risfu_mapping,
        reg_ids=reg_ids,
        nsoil=2,
    )

    HBdb = {
        1: DummyHB(subsurface=subsurface, cid_rank_mapping=cid_rank_mapping),
    }

    remote_payload = np.array([[0.7, 1.7]])
    print(f"[test-dz] preloading remote payload: {remote_payload}", flush=True)
    bus.preload(source=1, dest=0, tag=dz_tag_offset + 2001, payload=remote_payload)

    def fake_waitall(requests):
        print(f"[fake_waitall-dz] n_requests={len(requests)}", flush=True)
        for req in requests:
            req.complete()

    class _FakeMPI:
        class Request:
            @staticmethod
            def Waitall(requests):
                fake_waitall(requests)

    monkeypatch.setattr(mssubsurface, "MPI", _FakeMPI)

    exchange_dz_regional_units(cids=cids, rank=rank, HBdb=HBdb)
    print(f"[test-dz] sent_log entries={len(bus.sent_log)}", flush=True)

    sends = [
        entry
        for entry in bus.sent_log
        if entry[0] == 0 and entry[1] == 1 and entry[2] == dz_tag_offset + 1002
    ]
    print(f"[test-dz] matching sends count={len(sends)}", flush=True)
    assert len(sends) == 1
    print(f"[test-dz] send payload found={sends[0][3]}", flush=True)
    np.testing.assert_allclose(sends[0][3], np.array([[0.6, 1.6]]))

    expected_reg_dz = np.array(
        [
            [0.5, 1.5],
            [0.6, 1.6],
            [0.7, 1.7],
        ]
    )
    print(f"[test-dz] expected reg_dz_gw=\n{expected_reg_dz}", flush=True)
    print(f"[test-dz] actual reg_dz_gw=\n{HBdb[1].mssubsurface.reg_dz_gw}", flush=True)
    np.testing.assert_allclose(HBdb[1].mssubsurface.reg_dz_gw, expected_reg_dz)
    print("[test_exchange_dz_regional_units_external_exchange_and_assembly] done", flush=True)