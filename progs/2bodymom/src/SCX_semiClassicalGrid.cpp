
SCX_semiClassicalGrid::SCX_semiClassicalGrid(MeanFieldNucleusThick* nuc,Event &e,int bpoints,int zpoints) : 
	_nuc(nuc),
	_event(e),
	_bpoints(bpoints),
	_zpoints(zpoints),
	_bstep(0.),
	_zstep(0.),
	_pdens_fctr(nuc->getZ()),
	_ndens_fctr(nuc->getN()),
	_grid(NULL)
	{
	_bstep =    _nuc->getRange()/(_bpoints-1.); // 0 to nuc.getRange
	_zstep = 2.*_nuc->getRange()/(_zpoints-1.); // -nuc.getRange() to +nuc.getRange()
}
