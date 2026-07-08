#=

    Get metadata

=#


function get_Tg_offset(flight_number::Integer, ::Type{FT} = Float64, ::Val{use_summary_nc} = Val(false)) ::FT where {FT, use_summary_nc}
    if use_summary_nc
        summary = NC.Dataset(atlas_socrates_summary_file(flight_number), "r")
        flight_ind = findfirst(vec(Array(summary["flight_number"])) .== flight_number)
        return FT(summary[:deltaT][flight_ind])
    else

        if flight_number == 1
            return FT(1.44)
        elseif flight_number == 9
            return FT(-4.05)
        elseif flight_number == 10
            return FT(1.87)
        elseif flight_number == 11
            return FT(-1.48)
        elseif flight_number == 12
            return FT(-0.7)
        elseif flight_number == 13
            return FT(-0.95)
        else
            error("Flight number $flight_number not supported")
        end
    end
end



"""
    get_socrates_initial_time(flight_number)
Atlas LES simulation start time: `reference_time - 12h` from `SOCRATES_summary.nc` for `flight_number`.
"""
function get_socrates_initial_time(flight_number::Integer, ::Val{use_summary_nc} = Val(false)) where {use_summary_nc}
    if use_summary_nc
        summary = NC.Dataset(atlas_socrates_summary_file(flight_number), "r")
        flight_ind = findfirst(vec(Array(summary["flight_number"])) .== flight_number)
        return summary["reference_time"][flight_ind] - Dates.Hour(12)
    else
        if flight_number == 1
            return Dates.DateTime(2018, 1, 15, 14, 30, 0)
        elseif flight_number == 9
            return Dates.DateTime(2018, 2, 4, 15, 30, 0)
        elseif flight_number == 10
            return Dates.DateTime(2018, 2, 7, 13, 0, 0)
        elseif flight_number == 11
            return Dates.DateTime(2018, 2, 16, 16, 0, 0)
        elseif flight_number == 12
            return Dates.DateTime(2018, 2, 17, 16, 0, 0)
        elseif flight_number == 13
            return Dates.DateTime(2018, 2, 19, 15, 0, 0)
        else
            error("Flight number $flight_number not supported")
        end
    end
end