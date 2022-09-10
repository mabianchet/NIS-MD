import os
def Volfeat(irow):
    """Calculate the Volumes of the sites defined in tunnel.cfg as features"""
    import subprocess as sp
    nf,row=irow
    fn=row['traj_fn']
    fnxtc=tmp_dir+'/sel_'+os.path.basename(fn).replace('.dcd','.xtc')
    voldat=tmp_dir+'/volume_'+str(nf)+'.dat'
    tops=row['top_abs_fn']
    md.load(fn,top=tops).atom_slice(top_protein).save_xtc(fnxtc)
    process_string_list=['/opt/epock-1.0.5-Darwin-x86_64/bin/epock',
                         '-s','tmp.pdb',
                         '-c','Tunnels/tunnelf.cfg',
                         '-f',fnxtc,
                         '-o',voldat,
                         '--radii','Iodo.radii']
    try:
        sp.run(args=process_string_list,universal_newlines=True,
                       stdout=sp.PIPE, stderr=sp.PIPE)
        #os.remove(fnxtc)
        data= np.loadtxt(voldat, skiprows=1)
        try:
            data=np.delete(data,0,axis=1)
        except:
            data=np.delete(data,0)
        return nf,data
    except FileNotFoundError as e:
        print(e)
