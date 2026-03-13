using Test
import NCDatasets as NC
import SOCRATESSingleColumnForcings as SSCF

@testset "NCDatasets compatibility guards" begin
    @testset "Array(var) preserves shape" begin
        mktempdir() do dir
            path = joinpath(dir, "shape_test.nc")

            ds = NC.Dataset(path, "c")
            NC.defDim(ds, "x", 2)
            NC.defDim(ds, "y", 3)
            v = NC.defVar(ds, "v", Float64, ("x", "y"))
            vals = reshape(collect(1.0:6.0), 2, 3)
            v[:, :] = vals
            close(ds)

            dsr = NC.Dataset(path, "r")
            var = dsr["v"]
            @test length(var[:]) == 6
            @test size(Array(var)) == (2, 3)
            @test Array(var) == vals
            close(dsr)
        end
    end

    @testset "combine_air_and_ground_data handles insertion index shape" begin
        # vardata: [x, y, z, t]
        vardata = reshape(collect(1.0:24.0), 2, 1, 3, 4)
        # vardatag: [x, y, t] (no z dimension)
        vardatag = reshape(collect(101.0:108.0), 2, 1, 4)
        # insert index per [x, y, t], to be expanded inside combine_air_and_ground_data
        insert_location = fill(2, 2, 1, 4)

        out = SSCF.combine_air_and_ground_data(vardata, vardatag, 3; insert_location)

        @test size(out) == (2, 1, 4, 4)
        @test selectdim(out, 3, 2) == vardatag
    end
end
