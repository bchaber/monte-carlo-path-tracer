module Raytracer
    import Printf: @printf

    Tuple3 = Tuple{Float64, Float64, Float64}

    struct Hit
        l ::Float64
        d ::Tuple3
    end

    struct Scene
        c ::Tuple3
        r ::Float64
        ss ::Array{Scene, 1}
    end
    const sreps = sqrt(eps(1.0))
    corrected(x::Number) = floor(UInt8, (255.0clamp(x)^(1.0/2.2) + .5))
    normalize(v::Tuple3) = v./sqrt(dot(v,v))
    clamp(x::Number) = x < zero(x) ? zero(x) : x > one(x) ? one(x) : x
    dot(u::Tuple3,v::Tuple3) = sum(u .* v)

    function intersect(origin, direction, hit, sphere)
        v = sphere.c .- origin
        b = dot(v, direction)
        disc = b*b - dot(v, v) + sphere.r*sphere.r
        if disc < 0.
            disc = NaN
        else
            disc = sqrt(disc)
        end
        t1 = b - disc
        t2 = b + disc
        l  = (t2 > 0. ? (t1 > 0. ? t1 : t2) : Inf)
        if l > hit.l
            return hit
        else
            if size(sphere.ss, 1) == 0
                return Hit(l, normalize(origin .+ l .* direction .- sphere.c))
            else
                for scene in sphere.ss
                    hit = intersect(origin, direction, hit, scene)
                end
                return hit
            end
        end
    end

    function create(level, c, r)
        obj = Scene(c, r, [])
        if level == 1
            return obj
        else
            a = 3*r/sqrt(12)
            aux(x, z) = create(level-1, c .+ (x, a, z), r/2)
            Scene(c, 3*r, [obj, aux(-a, -a), aux(a, -a), aux(-a, a), aux(a, a)])
        end
    end

    const light = normalize((1.0, 3.0, -2.0))
    const ss = 2
    const n = 512
    const scene = create(2, (0.0, -1.0, 4.0), 1)
    const zero3 = (0.0, 0.0, 0.0)
    const hit0 = Hit(Inf, zero3)

    function raytrace(dir)
        hit = intersect(zero3, dir, hit0, scene)
        if hit.l == Inf
            return 0
        end
        g = dot(hit.d, light)
        if g < 0
            return 0.0
        end
        p = hit.l .* dir .+ sreps .* hit.d
        return intersect(p, light, hit0, scene).l < Inf ? 0.0 : g
    end

    aux(x, d) = x - n / 2.0 + d / ss
    
    open("image.pgm", "w") do f
      @printf(f, "P5 %d %d 255\n", n, n)
      for y in reverse(0:n-1)
        for x in reverse(0:n-1)
          g = 0.0
          for d in 0:ss^2-1
            g += raytrace(normalize((aux(x, d % ss), aux(y, d/ss), 1.0 * n)))
          end
          write(f, corrected(g/ss^2))
        end
      end
    end
end