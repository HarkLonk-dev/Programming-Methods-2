import turtle
import random
import math

SCREEN_WIDTH = 600
SCREEN_HEIGHT = 500

# ================= BASE =================
class GameObject(turtle.Turtle):
    def __init__(self, color, x, y, shape="square"):
        super().__init__()
        self.penup()
        self.shape(shape)
        self.color(color)
        self.goto(x, y)
        self.speed(0)
        self._start_x = x
        self._start_y = y

    def reset_pos(self):
        self.goto(self._start_x, self._start_y)
        self.setheading(0)
        self.showturtle()

# ================= OBJECTS =================
class Ball(GameObject):
    def __init__(self):
        super().__init__("white", 0, -140, shape="circle")
        self.shapesize(1.2)

class Keeper(GameObject):
    def __init__(self):
        super().__init__("#FFD700", 0, 150)
        self.shapesize(stretch_wid=2.5, stretch_len=1.8)

    def dive(self, target_x):
        error = random.randint(-40, 40)
        final_x = max(-180, min(180, target_x + error))
        self.goto(final_x, 150)

# ================= PLAYER =================
class PlayerVisual:
    def __init__(self):
        self.pen = turtle.Turtle()
        self.pen.hideturtle()
        self.pen.penup()
        self.reset()

    def draw_player(self, x, y):
        self.pen.clear()
        self.pen.goto(x, y)
        self.pen.color("#FF4500")
        self.pen.begin_fill()
        for _ in range(2):
            self.pen.forward(20); self.pen.left(90)
            self.pen.forward(35); self.pen.left(90)
        self.pen.end_fill()
        self.pen.goto(x + 10, y + 35)
        self.pen.dot(18, "#FFE4C4")

    def kick_anim(self):
        self.draw_player(-15, -140)

    def reset(self):
        self.draw_player(-60, -140)

# ================= UI =================
class TextManager:
    def __init__(self):
        self.score_pen = self._create_pen("white", -280, 220, ("Arial", 12, "bold"))
        self.center_pen = self._create_pen("yellow", 0, 40, ("Arial", 28, "bold"))
        self.help_pen = self._create_pen("white", 0, -220, ("Arial", 11, "italic"))

    def _create_pen(self, color, x, y, font_def):
        p = turtle.Turtle()
        p.hideturtle()
        p.penup()
        p.color(color)
        p.goto(x, y)
        p._font = font_def
        return p

    def update_ui(self, p1):
        self.score_pen.clear()
        self.score_pen.write(f"Score: {p1}", font=self.score_pen._font)

    def show_msg(self, text, color="yellow"):
        self.center_pen.clear()
        self.center_pen.color(color)
        self.center_pen.write(text, align="center", font=self.center_pen._font)

    def show_help(self, text):
        self.help_pen.clear()
        self.help_pen.write(text, align="center", font=self.help_pen._font)

# ================= GAME =================
class PenaltyGame:
    def __init__(self):
        self.screen = turtle.Screen()
        self.screen.bgcolor("#2E8B57")
        self.screen.setup(SCREEN_WIDTH, SCREEN_HEIGHT)
        self.screen.title("Penalty Game PRO")
        self.screen.tracer(0)

        self.draw_field()

        self.ui = TextManager()
        self.ball = Ball()
        self.keeper = Keeper()
        self.player = PlayerVisual()

        self.score = 0
        self.state = "PLAYING"

        # hướng mặc định
        self.dir_x = 0
        self.dir_y = 1

        # key
        self.screen.listen()
        self.screen.onkey(self.up, "Up")
        self.screen.onkey(self.left, "Left")
        self.screen.onkey(self.right, "Right")
        self.screen.onkey(self.kick, "space")

        self.ui.show_help("Arrow = aim | SPACE = shoot")
        self.screen.update()

    # ===== HƯỚNG =====
    def up(self):
        self.dir_x, self.dir_y = 0, 1
        self.ui.show_help("↑ Straight")

    def left(self):
        self.dir_x, self.dir_y = -4, 5
        self.ui.show_help("↖ Left")

    def right(self):
        self.dir_x, self.dir_y = 4, 5
        self.ui.show_help("↗ Right")

    # ===== SÂN =====
    def draw_field(self):
        p = turtle.Turtle()
        p.hideturtle()
        p.penup()
        p.color("white")
        p.pensize(3)
        p.goto(-180, 100)
        p.pendown()
        p.goto(-180, 180)
        p.goto(180, 180)
        p.goto(180, 100)
        p.penup()
        p.goto(0, -140)
        p.dot(8)

    # ===== SHOOT =====
    def kick(self):
        if self.state != "PLAYING":
            return

        self.state = "ANIMATING"
        self.player.kick_anim()

        target = self.dir_x * 80
        self.keeper.dive(target)

        self.animate_ball()

    # ===== ANIMATION (QUAN TRỌNG) =====
    def animate_ball(self):
        speed = 10    # 🔥 nhanh

        dx = self.dir_x
        dy = self.dir_y

        # chuẩn hóa vector
        length = (dx**2 + dy**2) ** 0.5
        dx /= length
        dy /= length

        if -300 < self.ball.xcor() < 300 and -200 < self.ball.ycor() < 200:
            self.ball.setx(self.ball.xcor() + dx * speed)
            self.ball.sety(self.ball.ycor() + dy * speed)

            self.screen.update()
            self.screen.ontimer(self.animate_ball, 20)
        else:
            self.check_result()

    # ===== RESULT =====
    def check_result(self):
        goal = abs(self.ball.xcor() - self.keeper.xcor()) > 30

        if goal:
            self.ui.show_msg("GOAL!!!", "lime")
            self.score += 1
        else:
            self.ui.show_msg("SAVED!", "red")

        self.ui.update_ui(self.score)

        self.screen.update()
        self.screen.ontimer(self.reset_turn, 1200)

    def reset_turn(self):
        self.ball.reset_pos()
        self.keeper.reset_pos()
        self.player.reset()
        self.state = "PLAYING"


        self.ui.show_help("Arrow = aim | SPACE = shoot")
        self.screen.update()

# RUN
PenaltyGame()
turtle.mainloop()