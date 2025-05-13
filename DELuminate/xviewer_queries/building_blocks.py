from typing import List

import pandas as pd

from configs.library_definition import BB_FILEPATH, BB_FILE_COLUMNS
from configs.xviewer import XViewer_Connection
from xviewer_queries.generic_query import generic_query


def query_candidate_bbs(
    bb_file_columns: List[str] = BB_FILE_COLUMNS,
    bb_permissions:str = 'X-Chem Pharmaceuticals',
    xviewer_connection:XViewer_Connection = None,
) -> pd.DataFrame:
    """
    Generates a dataframe of candidate BBs and filters the columns for only `bb_file_columns`
    """
    if xviewer_connection is not None:
        df = generic_query(
            query=f"""
            SELECT {','.join(bb_file_columns)} FROM XCRX_LIBREG.VW_XBBS
            WHERE
            ID IN
                (SELECT BB_ID
                FROM XCRX_LIBREG.LIBREG_BBPERMISSIONS
                WHERE
                PARTNER_ID IN (SELECT ID FROM XCRX_LIBREG.PRJM_PARTNERS WHERE PARTNER_NAME = '{bb_permissions}') AND
                NVL(END_DATE, SYSDATE+1) > SYSDATE)
            ORDER BY XBBID
            """,
            connection=xviewer_connection,
        )
    else:
        #Use a static export
        df = pd.read_csv(BB_FILEPATH).fillna("")
    df = df.loc[
        :,
        df.columns.intersection(bb_file_columns),
    ]
    return df

def lookup_integer_bbids(
    xbbids:List[str],
    xviewer_connection:XViewer_Connection,
    ) -> List[int]:
    """
    Converts a list of XBBIDs to a list of integer BBIDs
    """
    df = generic_query(
        """
        SELECT XBBID, ID as BBID 
        FROM XCRX_LIBREG.LIBREG_BBS 
        ORDER BY XBBID
        """,
        connection=xviewer_connection,
    )
    lookup_dict = pd.Series(df.BBID.values, index=df.XBBID).to_dict()
    try:
        bbids = [lookup_dict[xbbid] for xbbid in xbbids]
    except KeyError:
        failing_xbbids = []
        for xbbid in xbbids:
            try:
                _ = lookup_dict[xbbid]
            except KeyError:
                failing_xbbids.append(xbbid)
        raise Exception(f'The following XBBIDs did not have a corresponding BBID: {failing_xbbids}')
        
    return bbids

def lookup_external_bbids(
    xbbids:List[str],
    xviewer_connection:XViewer_Connection,
    ) -> List[str]:
    """
    Converts a list of XBBIDs to a list of external BBIDs ("EBBIDs")
    """
    df = generic_query(
        """
        SELECT XBBID, BB_KEY
        FROM XCRX_LIBREG.LIBREG_BBS_EXTERNAL_KEYS
        WHERE COMPANY = 'External IDs'
        ORDER BY XBBID
        """,
        connection=xviewer_connection,
    )
    lookup_dict = pd.Series(df.BB_KEY.values, index=df.XBBID).to_dict()
    ebbids = []
    for xbbid in xbbids:
        try:
            ebbids.append(lookup_dict[xbbid])
        except KeyError:
            ebbids.append("Unknown")
    
    return ebbids
