require 'matrix'
srand(0)

$timestep = 0.05
$i = 0
$counter = 0
class Sphere
	attr_accessor :p, :v, :r, :c
	def initialize(p, r, v, c)
		@p = p
		@v = v
		@r = r
		@c = c
	end

	def to_s
		f = 0.8
		retval = "material 0 0 0 #{@c[0]*f} #{@c[1]*f} #{@c[2]*f} #{@c[0]} #{@c[1]} #{@c[2]} 0 0 0 0 0 0 50 0 0\n" + "sphere #{$i} #{@p[0]} #{@p[1]} #{@p[2]} #{@r}"
		$i+=1
		return retval
	end

	def collides?(s)
		if (@p - s.p).r()**2 < (s.r + @r)**2
			return true
		else
			return false
		end
	end

	def tick
		@p += @v*$timestep
	end

	def mass
		@r**2
	end
end

def cross(a,b)
	x = a[1]*b[2] - a[2]*b[1]
	y = a[2]*b[0] - a[0]*b[2]
	z = a[0]*b[1] - a[1]*b[0]
	return Vector[x,y,z]
end

output = []
while output.length < 80
	p = Vector[rand()*2 - 1, rand()*2 - 1, rand()*2 - 1]
	c = Vector[rand(), rand(), rand()]
	v = Vector[rand()*2 - 1, rand()*2 - 1, rand()*2 - 1] * 0.2
	r = rand()*0.2
	try = Sphere.new(p, r, v, c)
	if output.find { |s| try.collides?(s) }
	else
		output << try
	end
end

while $counter < 1000
	scn = "camera 0 0 -6.5 0 0 1 0 1 0 50 0.1 100\n"
	output.each do |s|
		scn += s.to_s + "\n"
	end
	File.open("tmp.scn", "w") { |f| f.write(scn) }
	`./src/raypro tmp.scn art/output-#{sprintf("%04d", $counter)}.jpg -width 1024 -height 1024 -max_depth 4`
	puts $counter

	output.each do |s|
		s.tick
	end
	collided = {}
	output.each do |s|
		collideswith = output.find { |s2| s2.object_id != s.object_id and s.collides?(s2) }
		if collideswith and collided[s] != collideswith and collided[collideswith] != s
			collided[s] = collideswith
			puts "Collision"
			normalvector = (collideswith.p - s.p)
			normalvector = normalvector * (1.0/normalvector.r())

			temp = Vector[normalvector[0],normalvector[1],normalvector[2]+1]

			perpvector2 = cross(normalvector, temp)
			perpvector2 = perpvector2 * (1.0/perpvector2.r())
			perpvector1 = cross(normalvector, perpvector2)

			v1perp1 = s.v.inner_product(perpvector1)
			v1perp2 = s.v.inner_product(perpvector2)
			v2perp1 = collideswith.v.inner_product(perpvector1)
			v2perp2 = collideswith.v.inner_product(perpvector2)

			v1par = -s.v.inner_product(normalvector)
			v2par = -collideswith.v.inner_product(normalvector)

			s.v = perpvector1*v1perp1 + perpvector2*v1perp2 + normalvector*v1par
			collideswith.v = perpvector1*v2perp1 + perpvector2*v2perp2 + normalvector*v2par
		end

		if(s.p[0] < -1 || s.p[1] > 1)
			s.v = Vector[-s.v[0], s.v[1], s.v[2]]
		end
		if(s.p[1] < -1 || s.p[1] > 1)
			s.v = Vector[s.v[0], -s.v[1], s.v[2]]
		end
		if(s.p[2] < -1 || s.p[2] > 1)
			s.v = Vector[s.v[0], s.v[1], -s.v[2]]
		end

	end
	$counter += 1
	$i = 0
end



