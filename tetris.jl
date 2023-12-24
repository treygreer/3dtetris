using Printf
using LinearAlgebra
using Plots
using StaticArrays
#plotlyjs()

struct Piece
	name::String
	cells::Vector{Vector{Int64}}
	sticks::Vector{Vector{Int64}}
end

pieces = (                                                
	Piece("r2", [[1,1,1], [2,1,1], [3,1,1], [4,1,1], [3,2,1]],          [[1,2,3,4],[3,5]]),
	Piece("r1", [[3,1,1], [2,1,1], [1,1,1], [1,2,1], [1,3,1]],          [[1,2,3,4,5]]),
	Piece("r3", [[1,1,1], [2,1,1], [2,2,1], [1,2,1], [1,1,2]],          [[1,2,3,4,1,5]]),
	Piece("r4", [[2,1,1], [2,2,1], [1,2,1], [1,3,1], [1,3,2], [3,2,1]], [[1,2,3,4,5],[2,6]]),
	Piece("b1", [[1,1,1], [1,2,1], [2,2,1], [2,3,1], [2,3,2]],          [[1,2,3,4,5]]),
	Piece("b2", [[1,1,1], [2,1,1], [3,1,1], [3,2,1], [4,2,1]],          [[1,2,3,4,5]]),
	Piece("b3", [[1,1,1], [2,1,1], [3,1,1], [3,2,1], [3,3,1], [3,2,2]], [[1,2,3,4,5],[4,6]]),
	Piece("b4", [[3,1,2], [3,1,1], [2,1,1], [1,1,1], [1,2,1], [1,3,1]], [[1,2,3,4,5,6]]),
	Piece("y1", [[1,1,1], [2,1,1], [3,1,1], [3,2,1], [3,2,2]],          [[1,2,3,4,5]]),
	Piece("y2", [[1,1,1], [2,1,1], [3,1,1], [2,2,1], [2,2,2]],          [[1,2,3],[2,4,5]]),
	Piece("y3", [[3,1,1], [2,1,1], [2,2,1], [2,2,2], [1,2,2]],          [[1,2,3,4,5]]),
	Piece("y4", [[2,1,1], [2,2,1], [2,2,2], [2,3,2], [1,2,1], [3,2,1]], [[1,2,3,4],[5,2,6]])
)


const num_pieces = length(pieces)
print("$num_pieces pieces\n")

function xrot(θ)
	round.(Integer,
		[1 0      0
	     0 cos(θ) -sin(θ)
	     0 sin(θ)  cos(θ)])
end

function yrot(θ)
	round.(Integer, 
		[ cos(θ) 0 sin(θ)
	      0      1 0
	     -sin(θ) 0 cos(θ)])
end

function zrot(θ)
	round.(Integer, [cos(θ) -sin(θ) 0
	                 sin(θ)  cos(θ) 0
	                 0       0      1])
end

function legal_coords(c)
	(c[1]>=0 && c[1]<=3 &&
	 c[2]>=0 && c[2]<=3 &&
	 c[3]>=0 && c[3]<=3)
end

function bit(coords)
	if(legal_coords(coords))
		bit = coords[1] + 4*coords[2] + 16*coords[3]
		try
 			return UInt64(1) << bit
		catch e
			print("coords=", coords, ", bit=", bit, "\n")
			throw(DomainError("error in bit(coords)"))
		end
	else
		0
	end
end

function make_bits(p::Piece, rotation, offset)
	b = UInt64(0)
	for coords in p.cells
		b |= bit(rotation*coords+offset)
	end
	b
end

#print(make_bits(pieces[1], zrot(π), [3,3,-1]))

function print_bits(b::UInt64) 
	for z = 3:-1:0
		level = (b >> (16*z)) & 0xffff
		for y = 3:-1:0
			row = (level >> (4*y)) & 0xf
			for x = 0:3
				print((row>>x) & 0x1)
			end
			print(" ")
		end
		print("\n")
	end
end

rotations = []
for x_rotation in (xrot(0), xrot(π/2), xrot(π), xrot(3π/2))
	for other_rotation in (yrot(0), yrot(π/2), yrot(π), yrot(3π/2),
						   zrot(π/2), zrot(3π/2))
		rot = other_rotation * x_rotation
		@assert det(rot) == 1.0
		for prev_rot = rotations
			@assert prev_rot != rot
		end
		push!(rotations, rot)
	end
end
#print("found ", length(rotations), " unique rotations\n")
@assert length(rotations) == 24

function count_bits(bits::UInt64)
	count = 0
	for i=1:64
		if (bits & 1) == 1
			count += 1
		end
		bits >>= 1
	end
	return count 
end

struct Position
	bits::UInt64  # bit vector of positioned piece cells
	rotation::SMatrix{3,3,Int64}
	offset::SVector{3,Int64}
end

mutable struct SolutionPiece  
	piece_idx::Int64  # index into pieces
	num_cells::Int64  # number of cells in the piece
	positions_start_idx::Int64
	positions_ct::Int64
	solution_idx::Int64  # solution index
end

positions::Vector{Position} = []
solution::Vector{SolutionPiece} = []
position_idx = 0

for piece_idx in 1:num_pieces
	global position_idx
	piece = pieces[piece_idx]
	piece_length = length(piece.cells)
	solution_piece = SolutionPiece(piece_idx, length(piece.cells), position_idx, 0, 0)
	for rotation in (piece_idx==1 ? [rotations[1]] : rotations)
		for xoffset in -1:4
			for yoffset in -1:4
				for zoffset in -1:4
					offset = [xoffset, yoffset, zoffset]
					local bits
					try
						bits = make_bits(piece, rotation, offset)
					catch e
						print("whoops:\n")
						print("piece=", piece, "\n")
						print("rotation=", rotation, "\n")
						print("offset=", offset, "\n")
						throw(DomainError("error returned from make_bits"))
					end
					if count_bits(bits) == piece_length
						#print("  x=$xoffset y=$yoffset z=$zoffset: \n")
						#print_bits(bits)
						push!(positions, Position(bits, rotation, offset))
						solution_piece.positions_ct = solution_piece.positions_ct+1
						position_idx = position_idx+1
					end
				end
			end
		end
	end
	push!(solution, solution_piece)
end

sort!(solution, by = sp->sp.positions_ct)
for piece in solution
	print("$(pieces[piece.piece_idx].name): $(piece.positions_ct) possible positions\n");
end

function find_solution(so_far::UInt64, level::Int64)::Bool
	sp = solution[level]
	for solution_idx in 1:sp.positions_ct
		position_bits = positions[sp.positions_start_idx + solution_idx].bits
		if (so_far & position_bits) == 0
			solution[level].solution_idx = solution_idx
			new_so_far = so_far | position_bits
			if level == num_pieces
				@assert num_pieces < 12 || new_so_far == 0xffffffffffffffff
				return true
			else
				if find_solution(new_so_far, level+1) 
					return true
				end
			end
		end
	end
	return false
end

@time num_solutions = find_solution(UInt64(0), 1)
@time num_solutions = find_solution(UInt64(0), 1)

##### plot the solution
colors = Dict('r'=>"red", 'b'=>"blue", 'y'=>"yellow")
markers = [:x, :+, :o, :diamond]
println()
plot3d(background_color="grey", legend=false)
for sp in solution
	solution_idx = sp.solution_idx
	piece = pieces[sp.piece_idx]
	name = piece.name
	@printf "%s :  position %d of %d\n" name solution_idx sp.positions_ct
	position = positions[sp.positions_start_idx + solution_idx]
	rotation = position.rotation
	#offset = sp.positions[solution_idx].offset
	offset = position.offset
	color = colors[name[1]]
	marker = markers[Int(name[2])-Int('0')]
	for stick_element in piece.sticks
		x=[]; y=[]; z=[]
		for cell_idx in stick_element
			cell = piece.cells[cell_idx]
			coords = rotation*cell + offset
			push!(x, coords[1])
			push!(y, coords[2])
			push!(z, coords[3])
		end
		plot3d!(x,y,z, color=color, marker=marker, label=name)
	end
#		x = [(rotation*cell+offset)[1] for cell in piece.cells]
end
show(plot3d!())


##### try to find all of the solutions  (too slow)
println()
solution_ct::Int64 = 10
function count_solutions(so_far::UInt64, level::Int64)::Bool
	sp = solution[level]
	for position_idx in 1:sp.positions_ct
		position_bits = positions[sp.positions_start_idx + position_idx].bits
		if (so_far & position_bits) == 0
			new_so_far = so_far | position_bits
			if level == num_pieces
				#@assert num_pieces < 12 || new_so_far == 0xffffffffffffffff
				print(".")
				global solution_ct
				solution_ct = solution_ct - 1
				if solution_ct <= 0
					return true
				else
					return false
				end
			else
				if count_solutions(new_so_far, level+1) 
					return true
				end
			end
		end
	end
	return false
end
#solution_ct = 10
#@time count_solutions(UInt64(0), 1)
#solution_ct = 10
#@time count_solutions(UInt64(0), 1)
plot3d!()