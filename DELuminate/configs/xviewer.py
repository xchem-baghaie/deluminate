import getpass
import json

class XViewer_Connection():
    def __init__(self,
                 username:str = None,
                 password:str = None,
                 hostname:str = None,
                 port:int = None,
                 service_name:str = None,
                 ):
        with open('./DELuminate/configs/xviewer_connection.json') as f:
            cfg = json.load(f)
        self.username = username
        self.password = password
        if hostname is None:
            self.hostname = cfg['HOSTNAME']
        else:
            self.hostname = hostname
        if port is None:
            self.port = cfg['PORT']
        else:
            self.port = port
        if service_name is None:
            self.service_name = cfg['SERVICE_NAME']
        else:
            self.service_name = service_name
    
    def get_credentials(self) -> None:
        if self.username is None:
            self.username = getpass.getpass("Database username: ")
        if self.password is None:
            self.password = getpass.getpass("Database password: ")
        return None