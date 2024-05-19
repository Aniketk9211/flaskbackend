from server import create_app
import pytest

@pytest.fixture
def app():
    app = create_app()
    app.config.update({"TESTING": True})
    yield app

@pytest.fixture
def client(app):
    return app.test_client()

def test_convert_endpoint_with_valid_moldata(client):
    # Replace this with valid MOL data
    valid_moldata = '''CT1001419667


 12 12  0  1  0               999 V2000
   -0.0167    1.3781    0.0096 C   0 00  0  0  0  0  0  0  0  0  0  0
    0.0021   -0.0041    0.0020 C   0 00  0  0  0  0  0  0  0  0  0  0
    1.1709    2.0855    0.0021 C   0 00  0  0  0  0  0  0  0  0  0  0
    1.2084   -0.6789   -0.0132 C   0 00  0  0  0  0  0  0  0  0  0  0
    2.3960    0.0285   -0.0212 C   0 00  0  0  0  0  0  0  0  0  0  0
    2.3773    1.4107   -0.0131 C   0 00  0  0  0  0  0  0  0  0  0  0
   -0.9592    1.9054    0.0170 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9258   -0.5567    0.0083 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.1563    3.1654    0.0077 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.2231   -1.7588   -0.0184 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.3385   -0.4987   -0.0324 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.3051    1.9634   -0.0197 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2 02  0  1  0  0
  1  3 01  0  1  0  0
  1  7  1  0  0  0  0
  2  4 01  0  1  0  0
  2  8  1  0  0  0  0
  3  6 02  0  1  0  0
  3  9  1  0  0  0  0
  4  5 02  0  1  0  0
  4 10  1  0  0  0  0
  5  6 01  0  1  0  0
  5 11  1  0  0  0  0
  6 12  1  0  0  0  0
M  END
$$$$'''

    response = client.get(f"/convert?moldata={valid_moldata}")
    assert response.status_code == 200
    assert b"molecule.png" in response.data
    # Add more assertions as needed

def test_convert_endpoint_with_invalid_moldata(client):
    # Replace this with invalid MOL data
    invalid_moldata = '''


 12 12  0  1  0               999 V2000
   -0.0167    1.3781    0.0096 C   0 00  0  0  0  0  0  0  0  0  0  0
    0.0021   -0.0041    0.0020 C   0 00  0  0  0  0  0  0  0  0  0  0
    1.1709    2.0855    0.0021 C   0 00  0  0  0  0  0  0  0  0  0  0
    1.2084   -0.6789   -0.0132 C   0 00  0  0  0  0  0  0  0  0  0  0
    2.3960    0.0285   -0.0212 C   0 00  0  0  0  0  0  0  0  0  0  0
    2.3773    1.4107   -0.0131 C   0 00  0  0  0  0  0  0  0  0  0  0
   -0.9592    1.9054    0.0170 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9258   -0.5567    0.0083 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.1563    3.1654    0.0077 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.2231   -1.7588   -0.0184 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.3385   -0.4987   -0.0324 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.3051    1.9634   -0.0197 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2 02  0  1  0  0
  1  3 01  0  1  0  0
  1  7  1  0  0  0  0
  2  4 01  0  1  0  0
  2  8  1  0  0  0  0
  3  6 02  0  1  0  0
  3  9  1  0  0  0  0
  4  5 02  0  1  0  0
  4 10  1  0  0  0  0
  5  6 01  0  1  0  0
  5 11  1  0  0  0  0
  6 12  1  0  0  0  0
M  END
$$$$'''

    response = client.get(f"/convert?moldata={invalid_moldata}")
    assert response.status_code == 400
    assert b"Could not parse MOL data" in response.data
    # Add more assertions as needed
