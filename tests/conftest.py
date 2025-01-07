import boto3
import pytest
from moto import mock_aws


@pytest.fixture
def s3():
    with mock_aws():
        yield boto3.client("s3")
