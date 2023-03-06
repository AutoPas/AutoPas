function f()
    cutoff = 1.5
    li = [0.0, 0.0, 0.0]
    la = [3.0, 3.0, 3.0]
    dimensions = [3.0, 3.0, 3.0]
    for dimension in 2:3
        addArray = [0.0, 0.0, 0.0]
        addArray[dimension] = dimensions[dimension]
        upperMinCorner = [0.0, 0.0, 0.0]
        upperMaxCorner = [0.0, 0.0, 0.0]

        lowerMinCorner = [0.0, 0.0, 0.0]
        lowerMaxCorner = [0.0, 0.0, 0.0]
        for j in 1:3
            upperMaxCorner[j] = la[j] + cutoff
            lowerMinCorner[j] = li[j] - cutoff
            if j != dimension
                upperMinCorner[j] = li[j] - cutoff
                lowerMaxCorner[j] = la[j] + cutoff
            else
                upperMinCorner[j] = la[j]
                lowerMaxCorner[j] = li[j]
            end
        end
        println("corners for dimension: ", dimension)
        println("upper min: ", upperMinCorner, " | upper max: ", upperMaxCorner)
        println("lower min: ", lowerMinCorner, " | lower max: ", lowerMaxCorner)
    end
end

f()