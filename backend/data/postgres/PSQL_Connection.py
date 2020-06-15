import psycopg2
import PSQL_credentials as creds

class Postgres(object):
    """
    Estabilsh connection with a local chemblDB
    Requires a credential file with user password
    """
    _instance = None

    def __new__(cls):
        if cls._instance is None:
            cls._instance = object.__new__(cls)

            db_config = {'dbname': creds.PGDATABASE, 'host': creds.PGHOST,
                     'password': creds.PGPASSWORD, 'port': creds.pgPORT, 'user': creds.PGUSER}
            try:
                print('[Connecting to PostgreSQL database]')
                connection = Postgres._instance.connection = psycopg2.connect(**db_config)
                cursor = Postgres._instance.cursor = connection.cursor()
                cursor.execute('SELECT VERSION()')
                db_version = cursor.fetchone()

            except Exception as error:
                print('Error: connection not established {}'.format(error))
                Postgres._instance = None

            else:
                print('Connection established\n{}'.format(db_version[0]))

        return cls._instance

    def __init__(self):
        self.connection = self._instance.connection
        self.cursor = self._instance.cursor

    def query(self, query):
        try:
            result = self.cursor.execute(query)
        except Exception as error:
            print('error execting query "{}", error: {}'.format(query, error))
            return None
        else:
            return result

    def __del__(self):
        self.connection.close()
        self.cursor.close()


