class Sphere
	attr_accessor :x, :y, :z, :r
	def initialize(x, y, z, r)
		@x = x
		@y = y
		@z = z
		@r = r
	end

	def collides?(s)
		if (s.x - @x)**2 + (s.y - @y)**2 + (s.z - @z)**2 < (s.r + @r)**2
			return true
		else
			return false
		end
	end
end

output = []
while output.length < 300
	try = Sphere.new(rand()*2 - 1, rand()*2 - 1, rand()*2 - 1, (rand()*2 - 1)*0.2)
	if output.find { |s| try.collides?(s) }
	else
		output << try
	end
end

i = 0
output.each do |s|
	f = 0.75
	r = rand()
	g = rand()
	b = rand()
	puts "material 0 0 0 #{r*f} #{g*f} #{b*f} #{r} #{g} #{b} 0 0 0 0 0 0 65 0 0"
	puts "sphere #{i} #{s.x} #{s.y} #{s.z} #{s.r}"
	i+=1
end
