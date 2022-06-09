# -*- coding: utf-8 -*-
# snapshottest: v1 - https://goo.gl/zC4yUc
from __future__ import unicode_literals

from snapshottest import Snapshot
from snapshottest.file import FileSnapshot


snapshots = Snapshot()

snapshots['test_seqrun_stats 1'] = FileSnapshot('snap_test_seqrun_stats/test_seqrun_stats 1.tsv')
