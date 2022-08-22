use ggez::{event, graphics};
use ggez::{Context, ContextBuilder, GameResult};

use sepax2d::prelude::*;

mod dynamic;
use dynamic::DrawableShape;

const ANGLE: f32 = std::f32::consts::PI / 120.0;

trait DrawableRotation: DrawableShape + Rotate
{
}

impl DrawableRotation for Polygon
{
}

impl DrawableRotation for Capsule
{
}

impl DrawableRotation for Parallelogram
{
}

impl Shape for Box<dyn DrawableRotation>
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

struct MainState
{

    shapes: Vec<(Box<dyn DrawableRotation>, bool)>,

    mouse: (f32, f32),
    selected: Option<usize>

}

impl MainState
{

    fn new() -> GameResult<MainState>
    {

        let hexagon = Polygon::from_vertices((600.0, 100.0), vec![(0.0, 0.0), (50.0, 50.0), (50.0, 100.0), (0.0, 150.0), (-50.0, 100.0), (-50.0, 50.0)]);
        let triangle = Polygon::from_vertices((200.0, 100.0), vec![(0.0, 0.0), (100.0, 50.0), (0.0, 50.0)]);
    
        let vertical_capsule = Capsule::new((200.0, 400.0), (0.0, 30.0), 40.0);

        let gram = Parallelogram::new((600.0, 400.0), (25.0, 50.0), (25.0, -100.0));

        let s = MainState { shapes: vec!
        [

            (Box::new(hexagon), false),
            (Box::new(triangle), false),
            (Box::new(vertical_capsule), false),
            (Box::new(gram), false)

        ], mouse: (0.0, 0.0), selected: None };

        return Ok(s);

    }

}

impl event::EventHandler<ggez::GameError> for MainState
{

    fn update(&mut self, _context: &mut Context) -> GameResult
    {

        for (shape, collision) in self.shapes.iter_mut()
        {

            *collision = false;

            shape.rotate(ANGLE);

        }

        if let Some(i) = self.selected
        {

            self.shapes[i].0.set_position(self.mouse);

        }

        for i in 0..self.shapes.len()
        {

            for j in i+1..self.shapes.len()
            {

                if sat_overlap(&self.shapes[i].0, &self.shapes[j].0)
                {

                    self.shapes[i].1 = true;
                    self.shapes[j].1 = true;

                }

            }

        }

        return Ok(());

    }

    fn draw(&mut self, context: &mut Context) -> GameResult 
    {

        graphics::clear(context, [0.6, 0.6, 0.6, 1.0].into());

        for (shape, collides) in self.shapes.iter()
        {

            let mesh = shape.draw(*collides, context)?;
            let position = shape.position();

            graphics::draw(context, &mesh, ([position.0, position.1],))?;

        }

        let font = graphics::Font::new(context, "/PolandCanInto.otf")?;
        let text = graphics::Text::new(graphics::TextFragment
        {
            text: "Click and drag to move a shape! \n \n Shapes will change color when they are overlapping.".to_string(), 
            font: Some(font),
            color: Some(graphics::Color::new(0.0, 0.0, 0.0, 1.0)),
            ..Default::default()
        
        });

        graphics::draw(context, &text, ([25.0, 25.0],))?;

        graphics::present(context)?;

        return Ok(());

    }

    fn mouse_button_down_event(&mut self, _ctx: &mut Context, _button: ggez::input::mouse::MouseButton, x: f32, y: f32)
    {

        for (i, (shape, _overlap)) in self.shapes.iter().enumerate()
        {

            if contains_point(shape, (x, y))
            {

                self.selected = Some(i);
                return;

            }

        }

    }

    fn mouse_button_up_event(&mut self, _ctx: &mut Context, _button: ggez::input::mouse::MouseButton, _x: f32, _y: f32)
    {

        self.selected = None;

    }

    fn mouse_motion_event(&mut self, _ctx: &mut Context, x: f32, y: f32, _dx: f32, _dy: f32)
    {

        self.mouse = (x, y);

    }

}

fn main() -> GameResult
{

    let cb = ContextBuilder::new("rotation", "Meowitzher")
    .add_resource_path("./examples/resources/")
    .window_setup(ggez::conf::WindowSetup::default().title("Rotation Example!"));

    let (context, event_loop) = cb.build()?;
    let state = MainState::new()?;

    //Ommitted explicit return due to clippy complaining of unreachable code
    event::run(context, event_loop, state)

}