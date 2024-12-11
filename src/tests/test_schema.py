import unittest
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from src.config.database import Base
from src.database.models.schema import *  # Import all your models

class TestDatabaseSchema(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Create a test SQLite database
        cls.test_db_path = 'test_database.db'
        cls.engine = create_engine(f'sqlite:///{cls.test_db_path}')
        
        # Create all tables
        Base.metadata.create_all(cls.engine)
        
        # Create a session factory
        cls.Session = sessionmaker(bind=cls.engine)

    def setUp(self):
        # Create a new session for each test
        self.session = self.Session()

    def tearDown(self):
        # Close the session after each test
        self.session.close()

    @classmethod
    def tearDownClass(cls):
        # Remove the test database file
        if os.path.exists(cls.test_db_path):
            os.remove(cls.test_db_path)

    def test_create_and_query_basic_relationships(self):
        # Create test data
        accession = Accessions(
            ship_name="Test Ship",
            accession="TEST123",
            accession_tag="TEST",
            accession_new=1
        )
        self.session.add(accession)
        self.session.commit()

        ship = Ships(
            sequence="ATCG",
            md5="test_md5",
            accession_obj=accession
        )
        self.session.add(ship)
        self.session.commit()

        # Test relationships
        self.assertEqual(ship.accession_obj.ship_name, "Test Ship")
        self.assertEqual(accession.ships[0].sequence, "ATCG")

    # Add more specific tests for other relationships and constraints