#!/usr/bin/env python3
"""
Tests for the S3 utility functions.
"""
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch

import botocore.exceptions

from planter.database.utils.s3 import create_zip_archive, upload_to_s3


class TestS3Utils(unittest.TestCase):
    """Test cases for the S3 utility functions."""

    def setUp(self):
        """Set up test environment."""
        self.temp_dir = tempfile.mkdtemp()
        self.test_file_path = Path(self.temp_dir) / "test_file.txt"
        with open(self.test_file_path, "w") as f:
            f.write("Test content")

    def tearDown(self):
        """Clean up after tests."""
        import shutil

        shutil.rmtree(self.temp_dir)

    def test_create_zip_archive(self):
        """Test creating a zip archive."""
        zip_path = create_zip_archive(self.temp_dir)
        self.assertTrue(os.path.exists(zip_path))
        self.assertTrue(zip_path.endswith(".zip"))
        # Clean up the zip file
        os.remove(zip_path)

    @patch("boto3.client")
    def test_upload_to_s3_success(self, mock_boto3_client):
        """Test successful upload to S3."""
        # Mock the S3 client
        mock_s3 = MagicMock()
        mock_boto3_client.return_value = mock_s3

        # Create a proper ClientError for 404
        error_response = {"Error": {"Code": "404", "Message": "Not Found"}}
        mock_s3.head_object.side_effect = botocore.exceptions.ClientError(
            error_response, "HeadObject"
        )

        # Test the upload_to_s3 function
        result = upload_to_s3(self.temp_dir, "test_sample", "test_bucket")

        # Verify the result
        self.assertTrue(result)
        mock_s3.upload_file.assert_called_once()

    @patch("boto3.client")
    def test_upload_to_s3_file_exists(self, mock_boto3_client):
        """Test upload to S3 when file already exists."""
        # Mock the S3 client
        mock_s3 = MagicMock()
        mock_boto3_client.return_value = mock_s3

        # Mock the head_object method to not raise an exception (file exists)
        mock_s3.head_object.return_value = {}

        # Test the upload_to_s3 function
        result = upload_to_s3(self.temp_dir, "test_sample", "test_bucket")

        # Verify the result
        self.assertTrue(result)
        mock_s3.upload_file.assert_not_called()  # File should be skipped

    @patch("boto3.client")
    def test_upload_to_s3_error(self, mock_boto3_client):
        """Test upload to S3 with error."""
        # Mock the S3 client
        mock_s3 = MagicMock()
        mock_boto3_client.return_value = mock_s3

        # Create a proper ClientError for 404
        error_response = {"Error": {"Code": "404", "Message": "Not Found"}}
        mock_s3.head_object.side_effect = botocore.exceptions.ClientError(
            error_response, "HeadObject"
        )

        # Mock the upload_file method to raise an exception
        mock_s3.upload_file.side_effect = Exception("Upload error")

        # Test the upload_to_s3 function
        result = upload_to_s3(self.temp_dir, "test_sample", "test_bucket")

        # Verify the result
        self.assertFalse(result)
        mock_s3.upload_file.assert_called_once()


if __name__ == "__main__":
    unittest.main()
