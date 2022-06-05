use ggez::graphics::{Color, DrawMode, Mesh, MeshBuilder};
use ggez::{Context, GameResult};

use sepax2d::Shape;
use sepax2d::polygon::Polygon;
use sepax2d::circle::Circle;
use sepax2d::aabb::AABB;
use sepax2d::capsule::Capsule;

pub trait DrawableShape: Shape
{

    fn draw(&self, collides: bool, context: &mut Context) -> GameResult<Mesh>;

}

impl Shape for Box<dyn DrawableShape>
{

    fn position(&self) -> (f32, f32)
    {

        return self.as_ref().position();

    }

    fn set_position(&mut self, position: (f32, f32))
    {

        self.as_mut().set_position(position);

    }
    
    fn num_axes(&self) -> usize
    {

        return self.as_ref().num_axes();

    }

    fn get_axis(&self, index: usize, target: (f32, f32)) -> (f32, f32)
    {

        return self.as_ref().get_axis(index, target);

    }

    fn project(&self, axis: (f32, f32), normalize: bool) -> (f32, f32)
    {

        return self.as_ref().project(axis, normalize);

    }

    fn needs_closest(&self, index: usize) -> bool
    {

        return self.as_ref().needs_closest(index);

    }

    fn get_closest(&self, target: (f32, f32)) -> (f32, f32)
    {

        return self.as_ref().get_closest(target);

    }

    fn point(&self, index: usize) -> (f32, f32)
    {

        return self.as_ref().point(index);

    }

}

impl DrawableShape for Box<dyn DrawableShape>
{

    fn draw(&self, collides: bool, context: &mut Context) -> GameResult<Mesh>
    {

        return self.as_ref().draw(collides, context);

    }

}

impl DrawableShape for Polygon
{
    
    fn draw(&self, collides: bool, context: &mut Context) -> GameResult<Mesh>
    {

        //Very annoying that I don't know how to convert standard tuples
        //to GGEZ's glam tuples in a simple fashion, but I am too lazy to
        //find out how.
        let mut points = Vec::<[f32; 2]>::new();
        for (x, y) in self.vertices.iter()
        {

            points.push([*x, *y]);

        }

        let mesh = Mesh::new_polygon
        (

            context,
            DrawMode::fill(),
            &points,
            if collides { Color::BLUE } else { Color::WHITE }

        );

        return mesh;

    }

}

impl DrawableShape for Circle
{
    
    fn draw(&self, collides: bool, context: &mut Context) -> GameResult<Mesh>
    {
        
        let mesh = Mesh::new_circle
        (

            context,
            DrawMode::fill(),
            [0.0, 0.0],
            self.radius,
            1.0,
            if collides { Color::BLUE } else { Color::WHITE }

        );

        return mesh;

    }

}

impl DrawableShape for AABB
{
    
    fn draw(&self, collides: bool, context: &mut Context) -> GameResult<Mesh>
    {
        
        let mesh: GameResult<Mesh> = Mesh::new_rectangle
        (

            context,
            DrawMode::fill(),
            ggez::graphics::Rect { x: 0.0, y: 0.0, w: self.width, h: self.height },
            if collides { Color::BLUE } else { Color::WHITE }

        );

        return mesh;

    }

}

impl DrawableShape for Capsule
{
    
    fn draw(&self, collides: bool, context: &mut Context) -> GameResult<Mesh>
    {

        let arm = self.arm();
        let perp = self.perp();
        let radius = self.radius;

        let color = if collides { Color::BLUE } else { Color::WHITE };

        let mesh = MeshBuilder::new()
        .circle(DrawMode::fill(), [arm.0, arm.1], radius, 1.0, color)?
        .circle(DrawMode::fill(), [-arm.0, -arm.1], radius, 1.0, color)?
        .polygon(DrawMode::fill(), 
        &[

            [arm.0 + perp.0, arm.1 + perp.1],
            [arm.0 - perp.0, arm.1 - perp.1],
            [-arm.0 - perp.0, -arm.1 - perp.1],
            [-arm.0 + perp.0, -arm.1 + perp.1]

        ], color)?
        .build(context);

        return mesh;

    }

}