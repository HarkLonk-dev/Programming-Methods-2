import turtle
import random

# --- CẤU HÌNH HỆ THỐNG ---
SCREEN_WIDTH = 600
SCREEN_HEIGHT = 500

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

class Ball(GameObject):
    def __init__(self):
        super().__init__("white", 0, -140, shape="circle")
        self.shapesize(1.2)

class Keeper(GameObject):
    def __init__(self):
        super().__init__("#FFD700", 0, 150)
        self.shapesize(stretch_wid=2.5, stretch_len=1.8)

    def dive(self, target_x):
        error = random.randint(-70, 70) 
        final_x = max(-180, min(180, target_x + error))
        self.goto(final_x, 150)

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
        self.pen.color("black")
        self.pen.goto(x + 5, y); self.pen.pendown(); self.pen.goto(x + 5, y - 10); self.pen.penup()
        self.pen.goto(x + 15, y); self.pen.pendown(); self.pen.goto(x + 15, y - 10); self.pen.penup()

    def kick_anim(self):
        self.draw_player(-15, -140)

    def reset(self):
        self.draw_player(-60, -140)

class TextManager:
    def __init__(self):
        self.score_pen = self._create_pen("white", -280, 220, ("Arial", 12, "bold"))
        self.center_pen = self._create_pen("yellow", 0, 40, ("Arial", 28, "bold"))
        self.help_pen = self._create_pen("white", 0, -220, ("Arial", 11, "italic"))

    def _create_pen(self, color, x, y, font_def):
        p = turtle.Turtle(); p.hideturtle(); p.penup(); p.color(color); p.goto(x, y); p._font = font_def
        return p

    def update_ui(self, p1, p2, mode, turn):
        self.score_pen.clear()
        txt = f"MODE: {mode} | P1: {p1}" + (f" - P2: {p2} | Lượt: {turn}" if mode=="2P" else "")
        self.score_pen.write(txt, font=self.score_pen._font)

    def show_msg(self, text, color="yellow"):
        self.center_pen.clear(); self.center_pen.color(color)
        self.center_pen.write(text, align="center", font=self.center_pen._font)

    def show_help(self, text):
        self.help_pen.clear(); self.help_pen.write(text, align="center", font=self.help_pen._font)

class PenaltyGame:
    def __init__(self):
        self.screen = turtle.Screen()
        self.screen.bgcolor("#2E8B57")
        self.screen.setup(SCREEN_WIDTH, SCREEN_HEIGHT)
        self.screen.title("VGU Penalty Project - Team Flash")
        self.screen.tracer(0)

        self.draw_field()
        self.ui = TextManager()
        self.ball = Ball()
        self.keeper = Keeper()
        self.player = PlayerVisual()
        
        self.mode, self.p1_score, self.p2_score, self.round = None, 0, 0, 1
        self.is_p1, self.target_x, self.state = True, 0, "MENU"

        self.screen.listen()
        self.screen.onkey(self.sel_1p, "1"); self.screen.onkey(self.sel_2p, "2")
        self.screen.onkey(self.move_l, "Left"); self.screen.onkey(self.move_r, "Right")
        self.screen.onkey(self.kick, "space")
        
        self.ui.show_msg("PENALTY PRO", "white")
        self.ui.show_help("Nhấn '1' (1P) hoặc '2' (2P) để bắt đầu")
        self.screen.update()

    def draw_field(self):
        p = turtle.Turtle(); p.hideturtle(); p.penup(); p.color("white"); p.pensize(3)
        p.goto(-180, 100); p.pendown(); p.goto(-180, 180); p.goto(180, 180); p.goto(180, 100); p.penup()
        p.goto(0, -140); p.dot(8)

    def sel_1p(self): 
        if self.state == "MENU": self.mode="1P"; self.start()
    def sel_2p(self): 
        if self.state == "MENU": self.mode="2P"; self.start()

    def start(self):
        self.state = "PLAYING"; self.ui.center_pen.clear(); self.next_turn()

    def next_turn(self):
        limit = 5 if self.mode == "1P" else 10
        if self.round > limit: self.game_over(); return
        
        # Reset trạng thái cho lượt mới
        self.state = "PLAYING"
        self.ball.reset_pos()
        self.keeper.reset_pos()
        self.player.reset()
        self.ui.center_pen.clear() # Xóa chữ "VÀO" hoặc "HỤT" cũ
        
        self.target_x = 0
        self.draw_aim()
        self.ui.update_ui(self.p1_score, self.p2_score, self.mode, "P1" if self.is_p1 else "P2")
        self.ui.show_help("Mũi tên để ngắm - SPACE để sút")
        self.screen.update()

    def draw_aim(self):
        if not hasattr(self, 'aim'): self.aim = turtle.Turtle(); self.aim.hideturtle(); self.aim.penup()
        self.aim.clear(); self.aim.color("red"); self.aim.goto(self.target_x, 160); self.aim.dot(10); self.screen.update()

    def move_l(self):
        if self.state == "PLAYING" and self.target_x > -175: self.target_x -= 15; self.draw_aim()
    def move_r(self):
        if self.state == "PLAYING" and self.target_x < 175: self.target_x += 15; self.draw_aim()

    def kick(self):
        if self.state == "PLAYING":
            self.state = "ANIMATING"; self.aim.clear()
            self.player.kick_anim()
            self.keeper.dive(self.target_x)
            self.animate_ball()

    def animate_ball(self):
        if self.ball.ycor() < 150:
            self.ball.sety(self.ball.ycor() + 12)
            dx = (self.target_x - self.ball.xcor()) / 8
            self.ball.setx(self.ball.xcor() + dx)
            self.screen.update()
            self.screen.ontimer(self.animate_ball, 25)
        else:
            self.check_result()

    def check_result(self):
        goal = abs(self.ball.xcor() - self.keeper.xcor()) > 30
        if goal:
            self.ui.show_msg("VÀOOOOOO!", "#ADFF2F")
            if self.is_p1: self.p1_score += 1
            else: self.p2_score += 1
        else:
            self.ui.show_msg("BỊ CẢN PHÁ!", "#FF6347")
        
        if self.mode == "2P": self.is_p1 = not self.is_p1
        self.round += 1
        self.screen.update()
        # Chờ 1.5 giây rồi tự động gọi next_turn để reset vị trí
        self.screen.ontimer(self.next_turn, 1500)

    def game_over(self):
        self.state = "GAMEOVER"
        self.ui.show_msg("KẾT THÚC!", "white")
        self.screen.update()

if __name__ == "__main__":
    PenaltyGame()
    turtle.mainloop()