use ggez::{event, graphics};
use ggez::{Context, ContextBuilder, GameResult};

use sepax2d::prelude::*;

mod dynamic;
use dynamic::DrawableShape;

struct MainState {
    shapes: Vec<(Box<dyn DrawableShape>, bool)>,

    mouse: (f32, f32),
    selected: Option<usize>,
}

impl MainState {
    fn new() -> GameResult<MainState> {
        let hexagon = Polygon::from_vertices(
            (600.0, 100.0),
            vec![
                (0.0, 0.0),
                (50.0, 50.0),
                (50.0, 100.0),
                (0.0, 150.0),
                (-50.0, 100.0),
                (-50.0, 50.0),
            ],
        );
        let triangle =
            Polygon::from_vertices((150.0, 150.0), vec![(0.0, 0.0), (100.0, 50.0), (0.0, 50.0)]);

        let circle = Circle::new((300.0, 100.0), 80.0);
        let small_circle = Circle::new((700.0, 300.0), 40.0);

        let rectangle = AABB::new((550.0, 400.0), 100.0, 50.0);
        let square = AABB::new((100.0, 350.0), 50.0, 50.0);

        let vertical_capsule = Capsule::new((250.0, 400.0), (0.0, 30.0), 40.0);
        let rotated_capsule = Capsule::new((400.0, 500.0), (30.0, 20.0), 50.0);

        let gram = Parallelogram::new((400.0, 300.0), (25.0, 50.0), (25.0, -100.0));

        let s = MainState {
            shapes: vec![
                (Box::new(hexagon), false),
                (Box::new(triangle), false),
                (Box::new(circle), false),
                (Box::new(small_circle), false),
                (Box::new(rectangle), false),
                (Box::new(square), false),
                (Box::new(vertical_capsule), false),
                (Box::new(rotated_capsule), false),
                (Box::new(gram), false),
            ],
            mouse: (0.0, 0.0),
            selected: None,
        };

        return Ok(s);
    }
}

impl event::EventHandler<ggez::GameError> for MainState {
    fn update(&mut self, _context: &mut Context) -> GameResult {
        for (_shape, collision) in self.shapes.iter_mut() {
            *collision = false;
        }

        if let Some(i) = self.selected {
            self.shapes[i].0.set_position(self.mouse);
        }

        for i in 0..self.shapes.len() {
            for j in i + 1..self.shapes.len() {
                if sat_overlap(&self.shapes[i].0, &self.shapes[j].0) {
                    self.shapes[i].1 = true;
                    self.shapes[j].1 = true;
                }
            }
        }

        return Ok(());
    }

    fn draw(&mut self, context: &mut Context) -> GameResult {
        graphics::clear(context, [0.6, 0.6, 0.6, 1.0].into());

        for (shape, collides) in self.shapes.iter() {
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

    fn mouse_button_down_event(
        &mut self,
        _ctx: &mut Context,
        _button: ggez::input::mouse::MouseButton,
        x: f32,
        y: f32,
    ) {
        for (i, (shape, _overlap)) in self.shapes.iter().enumerate() {
            if contains_point(shape, (x, y)) {
                self.selected = Some(i);
                return;
            }
        }
    }

    fn mouse_button_up_event(
        &mut self,
        _ctx: &mut Context,
        _button: ggez::input::mouse::MouseButton,
        _x: f32,
        _y: f32,
    ) {
        self.selected = None;
    }

    fn mouse_motion_event(&mut self, _ctx: &mut Context, x: f32, y: f32, _dx: f32, _dy: f32) {
        self.mouse = (x, y);
    }
}

fn main() -> GameResult {
    let cb = ContextBuilder::new("overlap", "Meowitzher")
        .add_resource_path("./examples/resources/")
        .window_setup(ggez::conf::WindowSetup::default().title("Overlap Example!"));

    let (context, event_loop) = cb.build()?;
    let state = MainState::new()?;

    //Ommitted explicit return due to clippy complaining of unreachable code
    event::run(context, event_loop, state)
}
