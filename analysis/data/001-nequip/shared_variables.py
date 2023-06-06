import pathlib

path_to_traj = pathlib.Path('../../../runs/001-nequip/')    
path_to_output = {'rdf':pathlib.Path('rdf/'), 'msd':pathlib.Path('msd/'), 
    'elastic':pathlib.Path('elastic/'), 'pore':pathlib.Path('pore/'),
    'bad':pathlib.Path('bad/'),
    'bad_by_cn':pathlib.Path('bad_by_cn/'),
    'cn':pathlib.Path('cn/'),
    'reduction':pathlib.Path('reduced_trajectory/'),
    'ring':pathlib.Path('ring/')
    }
path_to_dataset = pathlib.Path("dataset/")
properties_dict = {
    'fast': {'rdf': True, 'msd': False,
                    'msd_by_run_id': False, 'elastic': True, 'pore': False, 'cn': True, 'bad': True,  'bad_by_cn':True,
                    'reduction': False, 'ring': False},
    'fewcores': {'rdf': True, 'msd': True,
                    'msd_by_run_id': True, 'elastic': True, 'pore': False, 'cn': False, 'bad': False,  'bad_by_cn':False,
                    'reduction': False, 'ring': False},                    
    'elastic': {'rdf': False, 'msd': False,
                        'msd_by_run_id': False, 'elastic': True, 'pore': False, 'cn': False, 'bad': False,  'bad_by_cn':False,
                        'reduction': False, 'ring': False},  
    'nopore': {'rdf': True, 'msd': True,
                    'msd_by_run_id': True, 'elastic': True, 'pore': False, 'cn': True, 'bad': True,  'bad_by_cn':True,
                    'reduction': True, 'ring': True}                         
}

runcategory = "001"