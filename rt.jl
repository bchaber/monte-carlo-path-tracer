module Raytracer
    import Printf: @printf

    Tuple1 = Tuple{Float64}
    Tuple3 = Tuple{Float64, Float64, Float64}

    struct Hit
        distance  :: Float64
        direction :: Tuple3
    end

    struct Sphere
        radius  :: Float64
        center  :: Tuple3
        emission :: Tuple1
        color :: Tuple1
    end
    corrected(x::Number) = floor(UInt8, (255.0clamp(x)^(1.0/2.2) + .5))
    normalize(v::Tuple3) = v./sqrt(dot(v,v))
    clamp(x::Number) = x < zero(x) ? zero(x) : x > one(x) ? one(x) : x
    cross(u::Tuple3, v::Tuple3) = (u[2]*v[3]-u[3]*v[2], u[3]*v[1]-u[1]*v[3],u[1]*v[2]-u[2]*v[1])
    dot(u::Tuple3,v::Tuple3) = sum(u .* v)

const lights = (
  Sphere(600., (50,681.6-.57,81.6),(0.25,), (.99,)),
  Sphere(3.,  (50,20.,40.),(0.25,), (.99,)),
)

const spheres = (
  Sphere(1e5,  (+1e5+01,40.8,81.6),(0.,), (.25,)),
  Sphere(1e5,  (-1e5+99,40.8,81.6),(0.,), (.25,)),
  Sphere(1e5,  (50,40.8, 1e5),     (0.,), (.99,)),
  Sphere(1e5,  (50,40.8,-1e5+170), (0.,), (.99,)),
  Sphere(1e5,  (50, 1e5, 81.6),    (0.,), (.99,)),
  Sphere(1e5,  (50,-1e5+81.6,81.6),(0.,), (.99,)),
  Sphere(16.5, (27,16.5,47),       (0.,), (.20,)),
  Sphere(16.5, (73,16.5,78),       (0.,), (.50,)),
)

    function intersect(origin, direction)
      t = 1e20
      hit = nothing
      for sph=(spheres..., lights...)
        d = intersect(origin, direction, sph)
        if d > 0. && d < t
          t = d
          hit = sph
        end
      end
      return t, hit
    end

    function intersect(origin, direction, sph::Sphere)
      op = sph.center .- origin
      b  = dot(op, direction)
      det = b*b - dot(op, op) + sph.radius*sph.radius
      if det < 0. return 0.0 end
      det = sqrt(det)
      if (b - det) > 1e-4 return b - det end
      if (b + det) > 1e-4 return b + det end
      return 0.
    end

    function raytrace(origin, direction)
      t, hit = intersect(origin, direction)
      if hit == nothing
          return (0.,)
      end
      f = (0.,)
      x = origin .+ t .* direction
      n = normalize(x .- hit.center)
      nl = dot(n, direction) < 0.0 ? n : -1.0.*n
      for lgh=lights
        vl = normalize(lgh.center .- x)
        t, hat = intersect(x, vl) 
        if hat == lgh
          f = f .+ hit.color .* hat.emission .* dot(nl, vl)
        end
      end
      return f
    end

    open("image.pgm", "w") do f
      w, h = 320, 240
      cam = (50.0, 50.0, 295.6)
      dir = (0.0, -0.042612, -1.0) |> normalize
      cx  = (0.5135w/h, 0.0, 0.0)
      cy  =  0.5135.*normalize(cross(cx, dir))
      @printf(f, "P5 %d %d 255\n", w, h)
      for y in reverse(1:h)
        for x in reverse(1:w)
          g = 0.0
          for sx=0:1
            for sy=0:1
              dx = 0.5
              dy = 0.5
              d  = @. cx.*( ( (sx + .5 + dx)/2.0 + x)/w - .5) .+
                      cy.*( ( (sy + .5 + dy)/2.0 + y)/h - .5) .+
                      dir;
              g += raytrace(cam .+ 140.0.*d, normalize(d))[1]
            end
          end
          write(f, corrected(g/4.0))
        end
      end
    end
end