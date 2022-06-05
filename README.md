# sepax2d
A safe Rust crate for finding and resolving collisions of convex polygons using the Separating Axis Theorem (SAT) in two dimensions.

### Usage

Add the following to the `[dependencies]` section of your `Cargo.toml`:

```toml
sepax2d = "0.3.0"
```

Import the types of shapes you'd like to use and create some new shapes:

```rust
use sepax2d::polygon::Polygon;
use sepax2d::circle::Circle;
use sepax2d::aabb::AABB; //Axis-aligned bounding boxes
use sepax2d::capsule::Capsule;

let rectangle = AABB::new(top_left, width, height);
let circle = Circle::new(center, radius);
let capsule = Capsule::new(center, arm, radius);

//The vertices of a polygon are position vectors relative to the shape's position, i.e.
//if position (1, 2), then the absolute location of the first vertex is (1, 0).
let triangle = Polygon::from_vertices(position, vec![(0.0, -2.0), (-1.0, 2.0), (1.0, 2.0)]);
```

Use the `sat_overlap(&left, &right)` method to get a `bool` indicating whether or not the two given shapes overlap.
Any struct implementing the `Shape` trait can be used.

```rust
use sepax2d::sat_overlap;

assert!(sat_overlap(&circle, &capsule));
assert!(sat_overlap(&triangle, &rectangle));
```

Use the `sat_collision(&left, &right)` method to get a `(f32, f32)` which represents the shift needed to add to `right`'s
position in order to resolve the overlap of the two shapes.

```rust
use sepax2d::{sat_overlap, sat_collision};

let resolution = sat_overlap(&circle, &triangle);

let position = triangle.position();
triangle.set_position((position.0 + resolution.0, position.1 + resolution.1));

assert!(!sat_overlap(&circle, &triangle));
```

Use the `contains_point(&shape, point)` method to get a `bool` indicating whether or not the specified point
is inside the given shape.

```rust
use sepax2d::aabb::AABB;
use sepax2d::contains_point;

let rectangle = AABB::new((0.0, 0.0), 5.0, 2.0);

assert!(contains_point(&rect, (1.0, 1.0)));
assert!(!contains_point(&rect, (10.0, -1.0)));
```

### Considerations
The Separating Axis Theorem only holds for convex shapes. If a concave shape is passed in, it is possible
for overlap to be missed and false overlap to be detected.

Convexity is not a problem for `Circle`, `AABB`, or `Capsule` as those are convex by definition. For 
polygons, it is possible to define a concave shape. Polygon convexity can be tested using the 
`polygon.is_convex()` method. Alternatively, the `Polygon::convex_from_vertices(position, vertices)`
constructor returns an `Option(Polygon)`, which is `None` if you try to make a concave polygon.

### Contribution
Please feel free to suggest additional features, bug fixes, or optimizations. Thanks!