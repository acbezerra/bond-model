module test

    function bond_pr_progress(v, vb, m, c, p,
                             r, gross_delta, iota,xi, k, alpha, pi, _lambda,
                             sigmal, sigmah, bondVmax; progress=nothing)
        if progress != nothing
            next!(progress)
        end
        return bond_pr_main(v, m, vb, vb, m, c, p,
                                 r, gross_delta, iota,xi, k, alpha, pi, _lambda,
                                 sigmal, sigmah, log(bondVmax/vb))
    end

    function pvvb_bond_surf(V0, pmin, pmax, pN, vN,
                            vbmin, vbmax, vbN,
                            m, c, p,
                            r, gross_delta, iota,
                            xi, k,
                            alpha, pi,
                            _lambda, sigmal, sigmah)

        bondVmax = get_bond_vmax(V0, m, c, p, sigmah, r,
                                 gross_delta, iota, xi, k,
                                 alpha, pi)

        pgrid = linspace(pmin, pmax, pN)
        vgrid = linspace(0, log(bondVmax/vbmin), vN)
        vbgrid = linspace(vbmin, vbmax, vbN)

        # cube_future = @spawn [bond_pr_main(v, m, vb, vb, m, c, p,
        #                          r, gross_delta, iota,xi, k, alpha, pi, _lambda,
        #                          sigmal, sigmah, log(bondVmax/vb)) for p=pgrid, v=vgrid, vb=vbgrid]

        p = Progress(length(pgrid) * length(vgrid) * length(vbgrid), 1)
        cube_future = @spawn [bond_pr_progress(v, vb, m, c, p,
                                 r, gross_delta, iota,xi, k, alpha, pi, _lambda,
                                 sigmal, sigmah, bondVmax; progress=p) for p=pgrid, v=vgrid, vb=vbgrid]
       return fetch(cube_future)

    end

end
