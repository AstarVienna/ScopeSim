import telescopy as tsp

def test_data_database():
    assert "Dave" in tsp.data.database.data_database()