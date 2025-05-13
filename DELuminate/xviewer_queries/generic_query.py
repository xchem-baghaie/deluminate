import oracledb
import pandas as pd

from configs.xviewer import XViewer_Connection

def generic_query(query:str, connection: XViewer_Connection)-> pd.DataFrame:
    with oracledb.connect(
        user = connection.username,
        password = connection.password,
        dsn=f'{connection.hostname}:{connection.port}/{connection.service_name}',
        disable_oob=True,
        ) as connection:
        with connection.cursor() as cursor:
            df = pd.DataFrame((cursor.execute(query)),columns=[x[0] for x in cursor.description])
    return df.fillna("")
