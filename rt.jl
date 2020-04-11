module Raytracer
    import Printf: @printf
    include("dual.jl")

    Float  = Float64#Dual{Float64}
    Tuple1 = Tuple{Float}
    Tuple3 = Tuple{Float, Float, Float}

    abstract type Material end
    struct Diffusive  <: Material end
    struct Specular   <: Material end

    struct Hit
        distance  :: Float
        direction :: Tuple3
    end

    struct Sphere
        radius  :: Float
        center  :: Tuple3
        emission :: Tuple1
        color :: Tuple1
        material :: Material
    end
    corrected(x::Dual) = 255.0clamp(x)^(1.0/2.2) + .5
    corrected(x::Number) = floor(UInt8, (255.0clamp(x)^(1.0/2.2) + .5))
    normalize(v) = v./sqrt(dot(v,v))
    clamp(x::Number) = x < zero(x) ? zero(x) : x > one(x) ? one(x) : x
    cross(u, v) = (u[2]*v[3]-u[3]*v[2], u[3]*v[1]-u[1]*v[3],u[1]*v[2]-u[2]*v[1])
    dot(u,v) = sum(u .* v)

const lights = (
  #Sphere(600., (50,681.6-.27,81.6),(0.25,), (.99,), Diffusive()),
  Sphere(1.,   (50.0,86.0,0.0),(0.25,), (.99,), Diffusive()),
)
const r = 1e5
const spheres = (
  Sphere(r,  (+r+01,   40.8,   81.6),(0.,), (.75,), Diffusive()),
  Sphere(r,  (-r+99,   40.8,   81.6),(0.,), (.25,), Diffusive()),
  Sphere(r,  (   50,   40.8,+r+0.00),(0.,), (.75,), Diffusive()),
  Sphere(r,  (   50,   40.8,-r+190.),(0.,), (.00,), Diffusive()),
  Sphere(r,  (   50,+r+0.00,   81.6),(0.,), (.75,), Diffusive()),
  Sphere(r,  (   50,-r+86.0,   81.6),(0.,), (.75,), Diffusive()),
  Sphere(9.5,(     9.5, 9.5, 9.5),(0.,), (.99,), Specular()),
  Sphere(9.5,(99.0-9.5, 9.5, 9.5),(0.,), (.99,), Diffusive()),
)

    function intersect(origin, direction)
      t = 1e20
      hit = nothing
      for sph=(spheres..., lights...)
        d = intersect(origin, direction, sph)
        if d > zero(d) && d < t*one(d)
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
      if det < zero(det) return zero(det) end
      det = sqrt(det)
      if (b - det) > 1e-4one(det) return b - det end
      if (b + det) > 1e-4one(det) return b + det end
      return zero(det)
    end

    function raytrace(origin, direction)
      t, hit = intersect(origin, direction)
      if hit ∈ lights
          return (one(t),)
      end
      if hit == nothing
          return (zero(t),)
      end
      f = (zero(t),)
      x = origin .+ t .* direction
      n = normalize(x .- hit.center)
      nl = dot(n, direction) < zero(t) ? n : -one(t).*n
      for lgh=lights
        vl = normalize(lgh.center .- x)
        t, hat = intersect(x, vl) 
        if hat == lgh
          f = f .+ hit.color .* hat.emission .* dot(nl, vl)
        end
      end
      return f
    end

    w, h = 267, 240
    c = zeros(w, h)
    ε = Dual(0.,1.)
    open("image.pgm", "w") do f
      cam = (50.0, 42.5, 295.6)
      dir = (0.0,   0.0, -1.0) |> normalize
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
              g += raytrace(cam .+ 145.0.*d, normalize(d))[1]
            end
          end
          #c[x, y] = corrected(g) |> partials
          write(f, corrected(.25g))
        end
      end
    end
end