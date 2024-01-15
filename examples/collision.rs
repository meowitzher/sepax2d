use ggez::{event, graphics};
use ggez::{Context, ContextBuilder, GameResult};

use sepax2d::aabb::AABB;
use sepax2d::capsule::Capsule;
use sepax2d::circle::Circle;
use sepax2d::polygon::Polygon;
use sepax2d::{contains_point, sat_collision, Shape};

mod dynamic;
use dynamic::DrawableShape;

struct MainState {
    shapes: Vec<(Box<dyn DrawableShape>, bool)>,

    mouse: (f32, f32),
    selected: Option<usize>,
}

const SPEED: f32 = 2.0;

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

        let circle = Circle::new((400.0, 200.0), 80.0);
        let small_circle = Circle::new((700.0, 300.0), 40.0);

        let rectangle = AABB::new((550.0, 400.0), 100.0, 50.0);
        let square = AABB::new((100.0, 350.0), 50.0, 50.0);

        let vertical_capsule = Capsule::new((250.0, 400.0), (0.0, 30.0), 40.0);
        let rotated_capsule = Capsule::new((400.0, 500.0), (30.0, 20.0), 50.0);

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
            ],
            mouse: (0.0, 0.0),
            selected: None,
        };

        return Ok(s);
    }
}

impl event::EventHandler<ggez::GameError> for MainState {
    fn update(&mut self, _context: &mut Context) -> GameResult {
        if let Some(i) = self.selected {
            let position = self.shapes[i].0.position();
            let new_position = (
                position.0 + SPEED * self.mouse.0,
                position.1 + SPEED * self.mouse.1,
            );

            self.shapes[i].0.set_position(new_position);

            let mut resolution = (0.0, 0.0);

            for (j, (shape, _collides)) in self.shapes.iter().enumerate() {
                if j != i {
                    let single_resolution = sat_collision(shape, &self.shapes[i].0);

                    resolution = (
                        resolution.0 + single_resolution.0,
                        resolution.1 + single_resolution.1,
                    );
                }
            }

            self.shapes[i]
                .0
                .set_position((new_position.0 + resolution.0, new_position.1 + resolution.1));
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
        let text = graphics::Text::new(graphics::TextFragment {
            text: "Click and hold to move a shape! \n \n Shapes will be stopped from overlapping."
                .to_string(),
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

    fn mouse_motion_event(&mut self, _ctx: &mut Context, _x: f32, _y: f32, dx: f32, dy: f32) {
        self.mouse = (dx, dy);
    }
}

fn main() -> GameResult {
    let cb = ContextBuilder::new("collision", "Meowitzher")
        .add_resource_path("./examples/resources/")
        .window_setup(ggez::conf::WindowSetup::default().title("Collision Example!"));

    let (context, event_loop) = cb.build()?;
    let state = MainState::new()?;

    //Ommitted explicit return due to clippy complaining of unreachable code
    event::run(context, event_loop, state)
}
