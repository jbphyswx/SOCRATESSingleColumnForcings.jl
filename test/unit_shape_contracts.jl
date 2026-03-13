using Test
import SOCRATESSingleColumnForcings as SSCF

@testset "Shape contract helpers" begin
    q_data = reshape(collect(1.0:12.0), 1, 1, 3, 4)

    @testset "read_profile_at_time preserves vertical shape" begin
        profile = SSCF.read_profile_at_time(q_data, 3, 4, 2)
        @test size(profile) == (3,)
        @test profile == [4.0, 5.0, 6.0]

        bad_q_data = reshape(collect(1.0:24.0), 2, 1, 3, 4)
        @test_throws Exception SSCF.read_profile_at_time(bad_q_data, 3, 4, 2)
    end

    @testset "read_profiles_over_time returns z-by-time" begin
        profiles = SSCF.read_profiles_over_time(q_data, 3, 4; time_indices = 2:4)
        @test size(profiles) == (3, 3)
        @test profiles == [4.0 7.0 10.0; 5.0 8.0 11.0; 6.0 9.0 12.0]

        bad_q_data = reshape(collect(1.0:24.0), 2, 1, 3, 4)
        @test_throws Exception SSCF.read_profiles_over_time(bad_q_data, 3, 4; time_indices = 1:4)
    end
end