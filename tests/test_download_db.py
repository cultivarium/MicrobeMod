import urllib.error
from unittest import mock

import pytest

from MicrobeMod import download_db


def http_error(code):
    return urllib.error.HTTPError(download_db.URL, code, "error", {}, None)


@pytest.fixture
def no_sleep(monkeypatch):
    sleep = mock.Mock()
    monkeypatch.setattr(download_db.time, "sleep", sleep)
    return sleep


@pytest.fixture(autouse=True)
def in_tmp_dir(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)


def test_client_error_fails_without_retry(no_sleep):
    with mock.patch.object(
        download_db.urllib.request, "urlopen", side_effect=http_error(403)
    ) as urlopen:
        with pytest.raises(SystemExit, match="HTTP 403"):
            download_db.download_database()
    assert urlopen.call_count == 1
    no_sleep.assert_not_called()


@pytest.mark.parametrize("code", [429, 503])
def test_retryable_http_error_is_retried(code, no_sleep, monkeypatch):
    monkeypatch.setattr(download_db, "MAX_ATTEMPTS", 3)
    with mock.patch.object(
        download_db.urllib.request, "urlopen", side_effect=http_error(code)
    ) as urlopen:
        with pytest.raises(SystemExit, match="after 3 attempts"):
            download_db.download_database()
    assert urlopen.call_count == 3
    assert no_sleep.call_count == 2


def test_sends_custom_user_agent(no_sleep):
    with mock.patch.object(
        download_db.urllib.request, "urlopen", side_effect=http_error(404)
    ) as urlopen:
        with pytest.raises(SystemExit):
            download_db.download_database()
    request = urlopen.call_args.args[0]
    assert request.get_header("User-agent") == download_db.USER_AGENT
