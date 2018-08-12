
function get_rgrow(svm)
   return rgrow(svm.pm.r, svm.pm.gross_delta) 
end


function get_cvm_vb(svm, sigma)
	return zhi_vb(svm.pm.m, svm.c, svm.p, 
                  sigma, svm.pm.r, svm.pm.gross_delta, 
                  svm.pm.iota, svm.pm.xi, svm.pm.kappa,
                  svm.pm.alpha, svm.pm.pi)
end


function get_param(svm, pname)
	val = NaN
	svm_pos = findin([string(x) for x in fieldnames(svm)], [pname])
	svm_pm_pos = findin([string(x) for x in fieldnames(svm.pm)], [pname])

	if pname == "C"
	    return svm.c * svm.pm.m
	elseif pname == "P"
	    return svm.p * svm.pm.m
	elseif pname == "delta"
	    return svm.pm.gross_delta - svm.pm.iota
	elseif length(svm_pos) > 0
	    return getfield(svm, fieldnames(svm)[svm_pos[1]])
	else
	    return getfield(svm.pm, fieldnames(svm.pm)[svm_pm_pos[1]])
	end
end


