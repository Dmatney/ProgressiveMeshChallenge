// Copyright 2016 Primal Space Systems, Inc. All rights reserved.

#pragma once

#pragma warning( disable : 4127 )		// conditional expression is constant
#pragma warning( disable : 4201 )		// nonstandard extension used: nameless struct/union
#pragma warning( disable : 4592 )		// symbol will be dynamically initialized (implementation limitation)

#include <cmath>
#include <cstdint>
#include <string>
#include <random>
#include <array>
#include <sstream>
#include <iomanip>
#include <assert.h>

#ifndef REAL_TYPE
#define REAL_TYPE double
#endif

/** Global fundamental types */
typedef int8_t int8;     // signed 8 bit integer aka char
typedef int16_t int16;   // signed 16 bit integer aka short
typedef int32_t int32;   // signed 32 bit integer
typedef int64_t int64;   // signed 64 bit integer

typedef uint8_t uint8;   // unsigned 8 bit integer aka byte
typedef uint16_t uint16; // unsigned 16 bit integer aka word
typedef uint32_t uint32; // unsigned 32 bit integer
typedef uint64_t uint64; // unsigned 64 bit integer

// bool                  // true or false
// float                 // 32 bit floating point number
// double                // 64 bit floating point number
// TCHAR                 // 2 byte UCS-2 character

// A small number for close to zero checks
#define COMMON_NEWLINE		"\n"
#define COMMON_EPSILON		0.000000001
#define EPSILON_BILLIONTH	0.0000000001
#define EPSILON_MILLIONITH	0.0000001
#define LOOSE_EPSILON		0.0001
#define PI					3.14159265359
#define SANE_VALUE			100000000.0

/** 
 * Compare two floating point numbers with an epsilon. 
 * Numbers greater than 1 are normalized by their magnitude
 */
inline bool EqualityComparison( const REAL_TYPE& Left, const REAL_TYPE& Right )
{
	if( fabs( Left ) < 1.0 )
	{
		return fabs( Left - Right ) < COMMON_EPSILON;
	}

	return fabs( ( Left - Right ) / Left ) < COMMON_EPSILON;
}

/**
 * Clamp a value (float, double, uint16, whatever) so that Min <= Value <= Max
 */
template <typename T>
T Clamp( const T Value, const T Min, const T Max )
{
	if( Value < Min )
	{
		return Min;
	}

	if( Value > Max )
	{
		return Max;
	}

	return Value;
}

/**
 * A vector of 3 x 32 bit signed integers. Used for specifying 3d coordinates in a grid or vertex indices of a triangle.
 */
class FIntVector
{
public:
	int32 X{ 0 };
	int32 Y{ 0 };
	int32 Z{ 0 };

	FIntVector() = default;

	FIntVector( int32 InX, int32 InY, int32 InZ )
		: X{ InX }
		, Y{ InY }
		, Z{ InZ }
	{
	}

	uint64 MemoryUsage() const
	{
		return sizeof( *this );
	}

	const std::string ToString() const
	{
		return "( " + std::to_string( X ) + ", " + std::to_string( Y ) + ", " + std::to_string( Z ) + " )";
	}

	bool IsValid() const
	{
		return X > -1 && Y > -1 && Z > -1;
	}

	// Make the values invalid so as to be sure to fire an assert if used
	void Invalidate()
	{
		X = -1;
		Y = -1;
		Z = -1;
	}

	int32& operator[] ( int32 Axis )
	{
		switch( Axis )
		{
		case 0:
			return X;
		case 1:
			return Y;
		case 2:
			return Z;
		default:
			assert( false );
			return X;
		}
	}

	int32 operator[] ( int32 Axis ) const
	{
		switch( Axis )
		{
		case 0:
			return X;
		case 1:
			return Y;
		case 2:
			return Z;
		default:
			assert( false );
			return X;
		}
	}

	bool operator<( const FIntVector& Other ) const
	{
		if( X == Other.X )
		{
			if( Y == Other.Y )
			{
				return( Z < Other.Z );
			}
			else
			{
				return( Y < Other.Y );
			}
		}

		return( X < Other.X );
	}

	inline int64 size() const
	{
		return 3;
	}
};

template <>
struct std::hash<FIntVector>
{
	std::size_t operator()( const FIntVector& Other ) const
	{
		return std::hash<int32>()( Other.X ) ^ std::hash<int32>()( Other.Y ) ^ std::hash<int32>()( Other.Z );
	}
};

/**
 * A vector of 3 floating point values. Used for specifying 3d point in space, a normal, a tangent, etc.
 */
union FVector
{
private:
	std::array<REAL_TYPE, 3> elem;

public:
	struct
	{
	REAL_TYPE X;
	REAL_TYPE Y;
	REAL_TYPE Z;
	};

	FVector()
		: X{ 0.0 }
		, Y{ 0.0 }
		, Z{ 0.0 }
	{
	}

	bool Validate() const
	{
		assert( fpclassify( X ) == FP_NORMAL || fpclassify( X ) == FP_ZERO );
		assert( fpclassify( Y ) == FP_NORMAL || fpclassify( Y ) == FP_ZERO );
		assert( fpclassify( Z ) == FP_NORMAL || fpclassify( Z ) == FP_ZERO );
		return true;
	}

	FVector( REAL_TYPE InX, REAL_TYPE InY, REAL_TYPE InZ )
		: X{ InX }
		, Y{ InY }
		, Z{ InZ }
	{
		assert( Validate() );
	}

	FVector( const FVector& Other )
		: X{ Other.X }
		, Y{ Other.Y }
		, Z{ Other.Z }
	{
	}

	FVector( FVector&& Other )
		: X{ Other.X }
		, Y{ Other.Y }
		, Z{ Other.Z }
	{
	}

	FVector& operator=( const FVector& Other )
	{
		X = Other.X;
		Y = Other.Y;
		Z = Other.Z;
		return *this;
	}

	uint64 MemoryUsage() const
	{
		return sizeof( *this );
	}

	const std::string ToString() const
	{
		return "( " + std::to_string( X ) + ", " + std::to_string( Y ) + ", " + std::to_string( Z ) + " )";
	}

	std::string ToHighPrecisionString() const
	{
		std::ostringstream stringStream;
		stringStream << std::setprecision(25);
		stringStream << "( " << X << ", " << Y << ", " << Z << " )";
		return stringStream.str();
	}

	std::string ToCPPCodeString() const
	{
		return "FVector" + ToHighPrecisionString();
	}

	// This operator is ONLY to be used for ordered set/map operations and has NO GEOMETRIC MEANING.
	bool operator<( const FVector& Other ) const
	{
		if( X < Other.X - COMMON_EPSILON )
		{
			return true;
		}
		else if( X > Other.X + COMMON_EPSILON )
		{
			return false;
		}

		if( Y < Other.Y - COMMON_EPSILON )
		{
			return true;
		}
		else if( Y > Other.Y + COMMON_EPSILON )
		{
			return false;
		}

		if( Z < Other.Z - COMMON_EPSILON )
		{
			return true;
		}
		else if( Z > Other.Z + COMMON_EPSILON )
		{
			return false;
		}

		return false;
	}

	REAL_TYPE& operator[] ( unsigned Axis )
	{
		assert(Axis < 3);
		return elem[Axis];
	}

	REAL_TYPE operator[] ( unsigned Axis ) const
	{
		assert(Axis < 3);
		return elem[Axis];
	}

	FVector operator +( const FVector& Other ) const
	{
		FVector Result
		{
			X + Other.X,
			Y + Other.Y,
			Z + Other.Z
		};

		assert( Result.Validate() );

		return Result;
	}

	FVector operator /(REAL_TYPE scalar) const
	{
		FVector Result
		{
			X / scalar,
			Y / scalar,
			Z / scalar
		};

		assert(Result.Validate());

		return Result;
	}

	FVector operator +=( const FVector& Other )
	{
		X += Other.X;
		Y += Other.Y;
		Z += Other.Z;

		assert( Validate() );
	
		return *this;
	}

	FVector operator -( const FVector& Other ) const
	{
		FVector Result
		{
			X - Other.X,
			Y - Other.Y,
			Z - Other.Z
		};

		assert( Result.Validate() );

		return Result;
	}

	FVector operator -=( const FVector& Other )
	{
		X -= Other.X;
		Y -= Other.Y;
		Z -= Other.Z;

		assert( Validate() );

		return *this;
	}

	FVector operator -() const
	{
		FVector Result { -X, -Y, -Z };

		assert( Result.Validate() );

		return Result;
	}

	FVector operator *( REAL_TYPE Other ) const
	{
		FVector Result
		{
			X * Other,
			Y * Other,
			Z * Other
		};

		assert( Result.Validate() );

		return Result;
	}

	FVector operator *=( REAL_TYPE Other )
	{
		X *= Other;
		Y *= Other;
		Z *= Other;

		assert( Validate() );

		return *this;
	}

	FVector operator *( const FVector& Other ) const
	{
		FVector Result
		{
			X * Other.X,
			Y * Other.Y,
			Z * Other.Z
		};

		assert( Result.Validate() );

		return Result;
	}

	FVector operator *=( const FVector& Other )
	{
		X *= Other.X;
		Y *= Other.Y;
		Z *= Other.Z;

		assert( Validate() );

		return *this;
	}

	bool operator ==(const FVector& Other) const
	{
		return ( EqualityComparison( X, Other.X )
				 && EqualityComparison( Y, Other.Y )
				 && EqualityComparison( Z, Other.Z ) );
	}

	bool operator !=(const FVector& Other) const
	{
		return ( !EqualityComparison( X, Other.X )
				 || !EqualityComparison( Y, Other.Y )
				 || !EqualityComparison( Z, Other.Z ) );
	}

	REAL_TYPE Dot( const FVector& Other ) const
	{
		return ( X * Other.X ) + ( Y * Other.Y ) + ( Z * Other.Z );
	}

	REAL_TYPE SquaredLength() const
	{
		return ( X * X ) + ( Y * Y ) + ( Z * Z );
	}

	REAL_TYPE Length() const
	{
		return sqrt( ( X * X ) + ( Y * Y ) + ( Z * Z ) );
	}

	REAL_TYPE Distance( const FVector& Other ) const
	{
		return ( *this - Other ).Length();
	}

	FVector Normalized() const
	{
		REAL_TYPE ReciprocalSize = ( REAL_TYPE )1.0 / Length();

		FVector Result
		{
			X * ReciprocalSize,
			Y * ReciprocalSize,
			Z * ReciprocalSize
		};

		assert( Result.Validate() );

		return Result;
	}

	FVector Cross( const FVector& Other ) const
	{
		FVector Result
		{
			Y * Other.Z - Z * Other.Y,
			Z * Other.X - X * Other.Z,
			X * Other.Y - Y * Other.X
		};

		assert( Result.Validate() );

		return Result;
	}

	FVector Abs() const
	{
		FVector Result
		{
			fabs( X ),
			fabs( Y ),
			fabs( Z )
		};

		assert( Result.Validate() );

		return Result;
	}

	int32 MaxComponent() const
	{
		FVector Other = Abs();

		if( Other.Z > Other.Y && Other.Z > Other.X )
		{
			return 2;
		}

		if( Other.Y > Other.X && Other.Y > Other.Z )
		{
			return 1;
		}

		return 0;
	}

	int32 MinComponent() const
	{
		FVector Other = Abs();

		if( Other.Z < Other.Y && Other.Z < Other.X )
		{
			return 2;
		}

		if( Other.Y < Other.X && Other.Y < Other.Z )
		{
			return 1;
		}

		return 0;
	}

	/**
	 * Returns the location of the point on the line
	 *
	 * @param [in] const FVector & StartLocation : one end of a line segment
	 * @param [in] const FVector & EndLocation : the other end of a line segment
	 *
	 * @retval bool : true if this point lies on the line segment.
	 */
	bool LocationOnLineSegment( const FVector& StartLocation, const FVector& EndLocation )
	{
		if( *this == StartLocation )
		{
			return true;
		}

		if( *this == EndLocation )
		{
			return true;
		}

		FVector LineSegment = EndLocation - StartLocation;
		FVector PointSegment = *this - StartLocation;

		REAL_TYPE CosineAngle = LineSegment.Normalized().Dot( PointSegment.Normalized() );
		if( fabs( CosineAngle - 1.0 ) > COMMON_EPSILON )
		{
			return false;
		}

		return LineSegment.Length() > PointSegment.Length();
	}

	static FVector randomUnit()
	{
		static std::random_device randomDevice;
		static std::default_random_engine randomEngine(randomDevice());
		static std::uniform_real_distribution<REAL_TYPE> uniformRandomnessSource(-1.0, 1.0);

		FVector randomVector;
		do 
		{
			// since we picked a random value along each axis this will be biased toward the corners of the unit cube
			// only accept vectors inside the unit sphere, so they will project uniformly on that sphere's surface
			randomVector = FVector( uniformRandomnessSource( randomEngine ), uniformRandomnessSource( randomEngine ), uniformRandomnessSource( randomEngine ) );
		}
		while( randomVector.SquaredLength() > 1.0 );

		return randomVector.Normalized();
	}

	FVector getOffsetInRandomDirection(REAL_TYPE inputOffsetDistance) const
	{
		//This function returns this given point offset by the given distance in a random direction
		return *this + randomUnit() * inputOffsetDistance;
	}

	FVector direction_to(const FVector& other_pt) const
	{
		// result: unit vector in the direction from this point to the other point passed in
		//         this *----->          * other

		FVector delta = other_pt - (*this);
		return delta.Normalized();
	}
};

template <>
struct std::hash<FVector>
{
	std::size_t operator()( const FVector& Other ) const
	{
		return std::hash<REAL_TYPE>()( Other.X ) ^ std::hash<REAL_TYPE>()( Other.Y ) ^ std::hash<REAL_TYPE>()( Other.Z );
	}
};

/** 
 * A class to represent a line segment
 */
class FLine
{
public:
	FVector Start;
	FVector End;
	bool bIsPoint{ true };

	FLine()
	{
	}

	FLine( const FVector& InStart, const FVector& InEnd )
		: Start{ InStart }
		, End{ InEnd }
	{
		bIsPoint = ( Start == End );
	}

	uint64 MemoryUsage() const
	{
		return sizeof( *this );
	}

	const std::string ToString() const
	{
		return Start.ToString() + " to " + End.ToString();
	}

	bool Validate() const
	{
		assert( Start.Validate() );
		assert( End.Validate() );
		return true;
	}

	FVector Delta() const
	{
		FVector Result
		{
			End.X - Start.X,
			End.Y - Start.Y,
			End.Z - Start.Z
		};

		assert( Result.Validate() );

		return Result;
	}
};

template <>
struct std::hash<FLine>
{
	std::size_t operator()( const FLine& Other ) const
	{
		return std::hash<FVector>()( Other.Start ) ^ std::hash<FVector>()( Other.End ) ^ std::hash<bool>()( Other.bIsPoint );
	}
};

/** 
 * A pair of vectors to represent an axis aligned bound box
 */
class FBoundingBox
{
public:
	FVector Mins;
	FVector Maxs;
	bool bInitialised;

	FBoundingBox()
		: Mins( REAL_TYPE( SANE_VALUE ), REAL_TYPE( SANE_VALUE ), REAL_TYPE( SANE_VALUE ) )
		, Maxs( REAL_TYPE( -SANE_VALUE ), REAL_TYPE( -SANE_VALUE ), REAL_TYPE( -SANE_VALUE ) )
		, bInitialised( false )
	{
	}

	FBoundingBox( const FVector& InMins, const FVector& InMaxs )
		: Mins{ InMins }
		, Maxs{ InMaxs }
		, bInitialised( true )
	{
	}

	uint64 MemoryUsage() const
	{
		return sizeof( *this );
	}

	const std::string ToString() const
	{
		return bInitialised ? ( Mins.ToString() + " -> " + Maxs.ToString() ) : "** Uninitialised **";
	}

	bool Validate() const
	{
		assert( Mins.Validate() );
		assert( Maxs.Validate() );

		assert( bInitialised );
		assert( Maxs.X >= Mins.X );
		assert( Maxs.Y >= Mins.Y );
		assert( Maxs.Z >= Mins.Z );
		return true;
	}

	const FVector Size() const
	{
		assert( Validate() );

		FVector Dimensions( Maxs - Mins );
		return Dimensions;
	}

	REAL_TYPE DiagonalLength() const
	{
		return Size().Length();
	}

	void AddPoint( const FVector& Point )
	{
		assert( Point.Validate() );

		// Check to see if we need to update the Mins
		if( Point.X < Mins.X )
		{
			Mins.X = Point.X;
		}

		if( Point.Y < Mins.Y )
		{
			Mins.Y = Point.Y;
		}

		if( Point.Z < Mins.Z )
		{
			Mins.Z = Point.Z;
		}

		// Check to see if we need to update the Maxs
		if( Point.X > Maxs.X )
		{
			Maxs.X = Point.X;
		}

		if( Point.Y > Maxs.Y )
		{
			Maxs.Y = Point.Y;
		}

		if( Point.Z > Maxs.Z )
		{
			Maxs.Z = Point.Z;
		}

		bInitialised = true;
	}

	void AddBounds( const FBoundingBox& Bounds )
	{
		assert( Bounds.Validate() );

		AddPoint( Bounds.Mins );
		AddPoint( Bounds.Maxs );
	}

	// Returns true if point touches or is inside bounding box
	bool Intersect( const FVector& inputPoint ) const
	{
		if( !bInitialised )
		{
			return false;
		}

		assert( Validate() );

		if( Maxs.X < inputPoint.X || Mins.X > inputPoint.X )
		{
			return false;
		}

		if( Maxs.Y < inputPoint.Y || Mins.Y > inputPoint.Y )
		{
			return false;
		}

		if( Maxs.Z < inputPoint.Z || Mins.Z > inputPoint.Z )
		{
			return false;
		}

		return true;
	}

	// Returns true if a pair of mins and maxs overlap in any way
	bool Intersect( const FBoundingBox& Other ) const
	{
		if( !bInitialised || !Other.bInitialised )
		{
			return false;
		}

		assert( Validate() );
		assert( Other.Validate() );

		if( Maxs.X < Other.Mins.X || Other.Maxs.X < Mins.X )
		{
			return false;
		}

		if( Maxs.Y < Other.Mins.Y || Other.Maxs.Y < Mins.Y )
		{
			return false;
		}

		if( Maxs.Z < Other.Mins.Z || Other.Maxs.Z < Mins.Z )
		{
			return false;
		}

		return true;
	}
};

template <>
struct std::hash<FBoundingBox>
{
	std::size_t operator()( const FBoundingBox& Other ) const
	{
		return std::hash<FVector>()( Other.Mins ) ^ std::hash<FVector>()( Other.Maxs ) ^ std::hash<bool>()( Other.bInitialised );
	}
};

/**
 * A vector of 4 floating point values. Used for specifying quaternions
 */
class FQuat
{
public:
	FVector Vector;
	REAL_TYPE W{ 1.0 }; // identity quaternion = 0,0,0,1

	FQuat() = default;

	FQuat( const FVector& InVec, REAL_TYPE InW )
		: Vector{ InVec }
		, W{ InW }
	{
	}

	FQuat( REAL_TYPE InX, REAL_TYPE InY, REAL_TYPE InZ, REAL_TYPE InW )
		: Vector{ InX, InY, InZ }
		, W{ InW }
	{
	}

	FQuat( const FVector& Angles )
	{
		REAL_TYPE CosHeading = cos( Angles.Y * ( REAL_TYPE )( PI / 360 ) );
		REAL_TYPE SinHeading = sin( Angles.Y * ( REAL_TYPE )( PI / 360 ) );
		REAL_TYPE CosAttitude = cos( Angles.Z * ( REAL_TYPE )( PI / 360 ) );
		REAL_TYPE SinAttitude = sin( Angles.Z * ( REAL_TYPE )( PI / 360 ) );
		REAL_TYPE CosBank = cos( Angles.X * ( REAL_TYPE )( PI / 360 ) );
		REAL_TYPE SinBank = sin( Angles.X * ( REAL_TYPE )( PI / 360 ) );

		REAL_TYPE CosHeadingCosAttitude = CosHeading * CosAttitude;
		REAL_TYPE SinHeadingSinAttitude = SinHeading * SinAttitude;

		Vector.X = ( CosHeadingCosAttitude * SinBank ) + ( SinHeadingSinAttitude * CosBank );
		Vector.Y = ( SinHeading * CosAttitude * CosBank ) + ( CosHeading * SinAttitude * SinBank );
		Vector.Z = ( CosHeading * SinAttitude * CosBank ) - ( SinHeading * CosAttitude * SinBank );
		W = ( CosHeadingCosAttitude * CosBank ) - ( SinHeadingSinAttitude * SinBank );
	}

	uint64 MemoryUsage() const
	{
		return sizeof( *this );
	}

	const std::string ToString() const
	{
		return "( ( " + std::to_string( Vector.X ) + ", " + std::to_string( Vector.Y ) + ", " + std::to_string( Vector.Z ) + " ), " + std::to_string( W ) + " )";
	}

	void Invert()
	{
		Vector = -Vector;
		// keep W the same
	}

	FQuat Inverse() const
	{
		return{ -Vector, W };
	}

	FVector operator*( const FVector& Point ) const
	{
		FVector Crossed = Vector.Cross( Point );
		return Point + ( Crossed * W * 2.0 ) + ( Vector.Cross( Crossed ) * 2.0 );
	}

	REAL_TYPE& operator[] ( unsigned Axis )
	{
		assert(Axis < 4);
		return (Axis < 3) ? Vector[Axis] : W;
	}

	REAL_TYPE operator[] ( unsigned Axis ) const
	{
		assert(Axis < 4);
		return (Axis < 3) ? Vector[Axis] : W;
	}
};

template <>
struct std::hash<FQuat>
{
	std::size_t operator()( const FQuat& Other ) const
	{
		return std::hash<FVector>()( Other.Vector ) ^ std::hash<REAL_TYPE>()( Other.W );
	}
};

/**
 * A class to contain the 3 coordinates of a triangle
 */
class FTriangle
{
public:
	std::array<FVector, 3> Corners;

	FTriangle() = default;

	FTriangle( const FVector& A, const FVector& B, const FVector& C )
	{
		Corners[0] = A;
		Corners[1] = B;
		Corners[2] = C;
	}

	uint64 MemoryUsage() const
	{
		return sizeof( *this );
	}

	bool IsDegenerate() const
	{
		FVector Edge1( Corners[1] - Corners[0] );
		FVector Edge2( Corners[2] - Corners[0] );

		if( Edge1.Cross( Edge2 ).Length() < COMMON_EPSILON )
		{
			return true;
		}

		return false;
	}

	FVector GetNormal() const
	{
		FVector Edge1(Corners[2] - Corners[0]);
		FVector Edge2(Corners[1] - Corners[0]);

		return Edge1.Cross( Edge2 );
	}

	const std::string ToString() const
	{
		return Corners[0].ToString() + " -> " + Corners[1].ToString() + " -> " + Corners[2].ToString();
	}

	FLine GetEdge( int32 Edge )
	{
		return FLine{ Corners[Edge % 3], Corners[( Edge + 1 ) % 3] };
	}
};

template <>
struct std::hash<FTriangle>
{
	std::size_t operator()( const FTriangle& Other ) const
	{
		return std::hash<FVector>()( Other.Corners[0] ) ^ std::hash<FVector>()( Other.Corners[1] ) ^ std::hash<FVector>()( Other.Corners[2] );
	}
};

/**
 * A class to contain and handle a plane equation (a normal and a signed distance)
 */
class FPlane
{
public:
	FVector Normal{ 0.0, 0.0, 1.0 }; // default = ground plane
	REAL_TYPE Distance{ 0.0 };

	FPlane() = default;

	bool operator ==( const FPlane& Other ) const
	{
		if( !EqualityComparison( Distance, Other.Distance ) )
		{
			return false;
		}
			
		return Normal == Other.Normal;
	}

	bool Validate() const
	{
		assert( Normal.Validate() );
		assert( fpclassify( Distance ) == FP_NORMAL || fpclassify( Distance ) == FP_ZERO );
		return true;
	}

	static FPlane FromTriangle( const FVector& Corner0, const FVector& Corner1, const FVector& Corner2 )
	{
		FPlane Plane;

		FVector Edge1( Corner2 - Corner0 );
		FVector Edge2( Corner1 - Corner0 );

		Plane.Normal = Edge1.Cross( Edge2 );
		Plane.Distance = -Plane.Normal.Dot( Corner0 );

		REAL_TYPE NormalMagnitude = Plane.Normal.Length();
		if( NormalMagnitude > COMMON_EPSILON )
		{
			const REAL_TYPE InverseMagnitude = ( ( REAL_TYPE ) 1.0 ) / NormalMagnitude;
			Plane.Normal *= InverseMagnitude;
			Plane.Distance *= InverseMagnitude;
		}

		assert( Plane.Validate() );
		return Plane;
	}

	static FPlane FromTriangle( const FTriangle& Triangle )
	{
		return FromTriangle( Triangle.Corners[0], Triangle.Corners[1], Triangle.Corners[2] );
	}

	FPlane( const FVector& NonUnitNormal, const FVector& Point )
		: Normal( NonUnitNormal )
		, Distance( -Normal.Dot( Point ) )
	{
		const REAL_TYPE inv_mag = ( ( REAL_TYPE ) 1.0 ) / Normal.Length();
		Normal *= inv_mag;
		Distance *= inv_mag;

		assert( Validate() );
	}

	// use values from already-computed plane coef's
	FPlane( REAL_TYPE A, REAL_TYPE B, REAL_TYPE C, REAL_TYPE D )
		: Normal { A, B, C }
		, Distance { D }
	{
		assert( Validate() );
	}

	FPlane operator -() const
	{
		FPlane Result{ -Normal.X, -Normal.Y, -Normal.Z, -Distance };

		assert( Result.Validate() );

		return Result;
	}

	uint64 MemoryUsage() const
	{
		return sizeof( *this );
	}

	const std::string ToString() const
	{
		return "( ( " + Normal.ToString() + " ), " + std::to_string( Distance ) + " )";
	}

	// SIGNED point-to-plane distance
	REAL_TYPE DistanceToPoint( const FVector& Point ) const
	{
		assert( Validate() );
		return Normal.Dot( Point ) + Distance;
	}

	REAL_TYPE DistanceToPointWithEpsilon( const FVector& Point ) const
	{
		REAL_TYPE Result = DistanceToPoint( Point );

		if( fabs( Result ) < COMMON_EPSILON )
		{
			// caller doesn't have to do its own epsilon check
			return 0.0;
		}

		return Result;
	}

	std::pair<int32, int32> ProjectedAxes() const
	{
		assert( Validate() );

		switch (Normal.MaxComponent())
		{
		case 0:
			return std::pair<int32, int32>(1, 2);
		case 1:
			return std::pair<int32, int32>(2, 0);
		case 2:
			return std::pair<int32, int32>(0, 1);
		default:
			assert( false );
			return std::pair<int32, int32>(0, 0);
		}
	}
};

template <>
struct std::hash<FPlane>
{
	std::size_t operator()( const FPlane& Other ) const
	{
		return std::hash<FVector>()( Other.Normal ) ^ std::hash<REAL_TYPE>()( Other.Distance );
	}
};

/** 
 * A class to manipulate a 3 by 3 transform matrix
 */
class FMatrix
{
public:
	std::array<FVector, 3> Elements;

	FMatrix Inverse() const
	{
		FMatrix Inverted;

		// computes the inverse of a matrix m
		REAL_TYPE Determinant = Elements[0][0] * (Elements[1][1] * Elements[2][2] - Elements[2][1] * Elements[1][2])
		                      - Elements[0][1] * (Elements[1][0] * Elements[2][2] - Elements[1][2] * Elements[2][0])
		                      + Elements[0][2] * (Elements[1][0] * Elements[2][1] - Elements[1][1] * Elements[2][0]);

		REAL_TYPE InverseDeterminant = ((REAL_TYPE) 1.0) / Determinant;

		Inverted.Elements[0][0] = (Elements[1][1] * Elements[2][2] - Elements[2][1] * Elements[1][2]) * InverseDeterminant;
		Inverted.Elements[0][1] = (Elements[0][2] * Elements[2][1] - Elements[0][1] * Elements[2][2]) * InverseDeterminant;
		Inverted.Elements[0][2] = (Elements[0][1] * Elements[1][2] - Elements[0][2] * Elements[1][1]) * InverseDeterminant;

		Inverted.Elements[1][0] = (Elements[1][2] * Elements[2][0] - Elements[1][0] * Elements[2][2]) * InverseDeterminant;
		Inverted.Elements[1][1] = (Elements[0][0] * Elements[2][2] - Elements[0][2] * Elements[2][0]) * InverseDeterminant;
		Inverted.Elements[1][2] = (Elements[1][0] * Elements[0][2] - Elements[0][0] * Elements[1][2]) * InverseDeterminant;

		Inverted.Elements[2][0] = (Elements[1][0] * Elements[2][1] - Elements[2][0] * Elements[1][1]) * InverseDeterminant;
		Inverted.Elements[2][1] = (Elements[2][0] * Elements[0][1] - Elements[0][0] * Elements[2][1]) * InverseDeterminant;
		Inverted.Elements[2][2] = (Elements[0][0] * Elements[1][1] - Elements[1][0] * Elements[0][1]) * InverseDeterminant;

		return Inverted;
	}

	uint64 MemoryUsage() const
	{
		return sizeof( *this );
	}
};

template <>
struct std::hash<FMatrix>
{
	std::size_t operator()( const FMatrix& Other ) const
	{
		return std::hash<FVector>()( Other.Elements[0] ) ^ std::hash<FVector>()( Other.Elements[1] ) ^ std::hash<FVector>()( Other.Elements[2] );
	}
};

