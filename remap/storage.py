import json
import sqlite3
from contextlib import contextmanager
from .errors import NotFoundError


SCHEMA = (
    "CREATE TABLE IF NOT EXISTS source_revisions("
    "id INTEGER PRIMARY KEY, name TEXT NOT NULL, payload TEXT NOT NULL)",
    "CREATE TABLE IF NOT EXISTS source_heads("
    "name TEXT PRIMARY KEY, source_id INTEGER NOT NULL REFERENCES source_revisions(id))",
    "CREATE TABLE IF NOT EXISTS target_revisions("
    "id INTEGER PRIMARY KEY, payload TEXT NOT NULL)",
    "CREATE TABLE IF NOT EXISTS plans("
    "id INTEGER PRIMARY KEY, name TEXT NOT NULL, "
    "source_id INTEGER NOT NULL REFERENCES source_revisions(id), "
    "target_id INTEGER NOT NULL REFERENCES target_revisions(id))",
    "CREATE TABLE IF NOT EXISTS results("
    "id INTEGER PRIMARY KEY, plan_id INTEGER NOT NULL REFERENCES plans(id), "
    "payload TEXT NOT NULL)",
    "CREATE TABLE IF NOT EXISTS receipts("
    "request_id TEXT PRIMARY KEY, plan_id INTEGER NOT NULL REFERENCES plans(id), "
    "result_id INTEGER NOT NULL UNIQUE REFERENCES results(id))",
)


def dumps(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False)


@contextmanager
def transaction(connection):
    connection.execute("BEGIN IMMEDIATE")
    try:
        yield
        connection.execute("COMMIT")
    except BaseException:
        connection.execute("ROLLBACK")
        raise


def connect(path):
    connection = sqlite3.connect(path, isolation_level=None)
    try:
        connection.execute("PRAGMA foreign_keys=ON")
        connection.execute("PRAGMA busy_timeout=5000")
        with transaction(connection):
            for sql in SCHEMA:
                connection.execute(sql)
        return connection
    except BaseException:
        connection.close()
        raise


def fetch_one(connection, sql, args, label):
    row = connection.execute(sql, args).fetchone()
    if row is None:
        raise NotFoundError(label)
    return row
