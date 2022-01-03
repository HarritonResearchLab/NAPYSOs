import pandas as pd
# import numpy as np
import sys


def main(i_f: str, m_f: str, o_f: str) -> pd.DataFrame:
    infrareds = pd.read_csv(i_f)
    updated_classes = pd.read_csv(m_f)

    # infrareds.set_index('ID_1')
    # infrareds.sort_index(inplace=True)
    # updated_classes.set_index('ID')
    # updated_classes.sort_index(inplace=True)

    u_oids = list(updated_classes['ID'])
    u_pclasses = list(updated_classes['primary_class'])

    o_oids = list(infrareds['ID_1'])
    o_pclasses = list(infrareds['primary_class'])

    for (u_oid, u_pclass) in zip(u_oids, u_pclasses):
        exists = False
        for (o_oid, o_pclass) in zip(o_oids, o_pclasses):
            if u_oid == o_oid:
                infrareds[infrareds['ID_1'] == o_oid]['primary_class'] = u_pclass
                exists = True
                break
        if not exists:
            pass  # handle this case if needed


    # infrareds[infrareds['ID_1'] == updated_classes['ID']]['primary_class'] = updated_classes[infrareds['ID_1'] == updated_classes['ID']]['primary_class']

    infrareds.to_csv(o_f)

    return infrareds


if __name__ == "__main__":
    print("outdated function--uncomment 'main' to run", file=sys.stderr)
    # main(i_f='objects_infrared.csv', m_f='merged_results.csv', o_f='updated_object_infrared.csv')
