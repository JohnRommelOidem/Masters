from manim_cnf import *
from manim import config

config.assets_dir = "assets"
default_camera = {
    'phi': 75*DEGREES,
    'theta': 30*DEGREES,
    "focal_distance": 1000
}


class TitleScene(Scene):
    def construct(self):
        title = VGroup(
            Text("Chaos and Escapes of Charged Particles", t2c={'[21:38]': BLUE}),
            Text("Orbiting Uniformly Magnetized", t2c={'[9:29]': YELLOW}),
            Text("Gravitational Centers").set_color(ORANGE)
        ).arrange(DOWN)
        by = VGroup(
            Text('By: John Rommel Oidem'),
            Text('Adviser: Michael Francis Ian Vega II')
        ).scale(0.5).arrange(DOWN).next_to(title, DOWN)
        up_logo = ImageMobject("UP_logo").scale(0.5).next_to(title, 1.2*UP)
        gravity_logo = ImageMobject("gravity_logo.png").scale(0.3).next_to(by, DOWN).shift(RIGHT)
        ganap_logo = ImageMobject("ganap_logo.png").scale(0.13).next_to(gravity_logo, LEFT)
        self.play(Write(title), Write(by), FadeIn(up_logo), FadeIn(gravity_logo), FadeIn(ganap_logo), run_time=0.5)
        self.wait(0.2)
        self.next_section()
        self.play(Unwrite(title), Unwrite(by), FadeOut(up_logo), FadeOut(gravity_logo), FadeOut(ganap_logo), run_time=0.5)
        self.wait(0.2)


class AstroScene(Scene):
    def construct(self):
        M87_image = ImageMobject("M87_image.png").scale(2)
        M87_image.add(Text("Messier 87").next_to(M87_image, DOWN))
        M87_polar = Group(ImageMobject("M87_polar.png").scale(2))
        M87_polar.add(Text("Messier 87").next_to(M87_polar, DOWN))
        SA_image = ImageMobject("SA_image.png").scale(2).to_edge(RIGHT)
        SA_image.add(Text("Sagittarius A*").next_to(SA_image, DOWN))
        SA_polar = ImageMobject("SA_polar.png").scale(2).to_edge(RIGHT)
        SA_polar.add(Text("Sagittarius A*").next_to(SA_polar, DOWN))
        self.play(FadeIn(M87_image))
        self.wait(0.2)
        self.next_section()
        self.play(FadeIn(M87_polar))
        self.remove(M87_image)
        self.wait(0.2)
        self.next_section()
        self.play(M87_polar.animate.to_edge(LEFT))
        self.play(FadeIn(SA_image))
        self.wait(0.2)
        self.next_section()
        self.play(FadeIn(SA_polar))
        self.remove(SA_image)
        self.wait(0.2)
        self.next_section()
        self.play(FadeOut(SA_polar, M87_polar[1]))
        jets = ImageMobject("jet.png")
        M87_polar.set_z_index(jets.z_index+1)
        text = Text('Astrophysical jets').to_corner(DL)
        self.play(
            FadeIn(jets), Write(text, run_time=0.5),
            FadeOut(M87_polar[0], scale=0.1, target_position=0.55 * UP + 0.25 * LEFT)
        )
        self.wait(0.2)
        self.next_section()
        line = Line(6*UP, 6*DOWN).rotate(32.1*DEGREES).shift(1.3*UP+0.6*LEFT).set_color(BLUE)
        line_out = Line(line.get_end(), line.get_start()).set_color(line.get_color())
        jet_list = BulletedList(
            'Linear structures',
            'Accretion disks',
            'Plasma outflow'
        ).to_corner(UR)
        self.play(Create(line, run_time=0.5))
        self.add(line_out)
        self.remove(line)
        self.play(Uncreate(line_out), Write(jet_list[0]), run_time=0.5)
        self.wait(0.2)
        self.next_section()
        disk = Ellipse(width=3.5, height=0.3).shift(0.55*UP+0.25*LEFT).rotate(31.5*DEGREES).set_color(BLUE)
        disk_out = Ellipse(width=-3.5, height=0.3).shift(0.55*UP+0.25*LEFT).rotate(31.5*DEGREES).set_color(disk.get_color())
        self.play(Create(disk))
        self.add(disk_out)
        self.remove(disk)
        self.play(Uncreate(disk_out), Write(jet_list[1], run_time=0.5))
        self.wait(0.2)
        self.next_section()
        up_arrow = Arrow(UP, 3*UP)
        down_arrow = Arrow(DOWN, 3*DOWN)

        line = VGroup(up_arrow, down_arrow).shift(0.55*UP+0.15*LEFT).rotate(32.1*DEGREES).set_color(BLUE)
        self.play(FadeIn(line[0], shift=PolarPlane().polar_to_point(2, (32.1+90)*DEGREES)),
                  FadeIn(line[1], shift=PolarPlane().polar_to_point(2, (32.1-90)*DEGREES)))
        self.play(Write(jet_list[2], run_time=0.5))
        self.wait(0.2)
        self.next_section()
        self.play(FadeOut(line[0], shift=PolarPlane().polar_to_point(3, (32.1+90)*DEGREES)),
                  FadeOut(line[1], shift=PolarPlane().polar_to_point(3, (32.1-90)*DEGREES)))
        self.wait(0.2)


class BHScene(ThreeDScene):
    def construct(self):
        jets = ImageMobject(f'jet')
        text = Text('Astrophysical jets').to_corner(DL)
        BH = BlackHole()
        jet_list = BulletedList(
            'Linear structures',
            'Accretion disks',
            'Plasma outflow'
        ).to_corner(UR)
        group = Group(jets, text, jet_list)
        self.set_camera_orientation(phi=90*DEGREES, theta=0*DEGREES)
        self.add_fixed_in_frame_mobjects(group, BH)
        self.remove(BH)
        self.wait(0.1)
        self.play(group.animate.scale(2).rotate((-32.8) * DEGREES).shift(1.3 * DOWN + 0.3 * LEFT))
        self.play(FadeOut(group, scale=3, shift=2 * DOWN + 0.5 * LEFT),
                  GrowFromCenter(BH))
        self.wait(0.1)
        self.next_section()
        field = ArrowVectorField(lambda _: np.array([0., 0., 2.]),
                                 x_range=[-7, 7, 2],
                                 y_range=[-7, 7, 2],
                                 z_range=[-7, 7, 2],
                                 color=YELLOW)
        self.play(Create(field))
        self.wait(0.2)
        self.next_section()
        self.play(Uncreate(field))
        L = -3.
        r_eq = np.roots([1, 0, 0, 4, -4 * L ** 2])
        r0 = np.abs(r_eq[np.logical_and(r_eq >= 0, np.isreal(r_eq))])[0]+0.3
        test = MNCTestParticle3D(L=L, E=1.2, r=r0, z=0, theta=np.pi, pr=-0.4)
        self.remove(BH)
        test.move()
        self.add(test, BH)
        self.wait(10)


class SystemIntroScene(ThreeDScene):
    def construct(self):
        BH = BlackHole()
        self.add(BH)
        BH_group = VGroup(
            Text("Schwarzschild").set_color(ORANGE),
            Text("black hole").set_color(ORANGE),
            BlackHole(),
            Text("Magnetized"),
            Text("black hole (MBH)")
        ).arrange(DOWN).to_edge(LEFT, MED_LARGE_BUFF)
        MNC_group = VGroup(
            Text("Point mass w/").set_color(ORANGE),
            Text("Newtonian gravity").set_color(ORANGE),
            MNCenter().rotate(-PI/2, RIGHT),
            Text("Magnetized Newtonian"),
            Text("center (MNC)")
        ).arrange(DOWN).to_edge(RIGHT, MED_LARGE_BUFF)
        self.play(Write(BH_group[:2], run_time=0.5), ReplacementTransform(BH, BH_group[2]))
        self.wait(0.2)
        self.next_section()
        self.play(Write(BH_group[3:], run_time=0.5))
        self.wait(0.2)
        self.next_section()
        self.play(Write(MNC_group[:2], run_time=0.5), GrowFromCenter(MNC_group[2]))
        self.wait(0.2)
        self.next_section()
        self.play(Write(MNC_group[3:], run_time=0.5))
        self.wait(0.2)
        self.next_section()
        self.play(FadeOut(BH_group), FadeOut(MNC_group))
        self.wait(0.2)
        self.next_section()


class ProblemStatement(Scene):
    def construct(self):
        problem_statemetn = BulletedList(
            "MBH has been studied aplenty",
            "We have found no study on the MNC",
            "Start with the simpler model, then build from there",
            "We focus on two aspects of both systems",
            "Escape dynamics",
            "Chaotic dynamics"
        )
        problem_statemetn[-2:].shift(RIGHT)
        for i in problem_statemetn:
            self.play(FadeIn(i, shift=UP))
            self.wait(0.2)
            self.next_section()
        self.play(FadeOut(problem_statemetn, shift=UP))
        question = BulletedList(
            r"What determines if a particle\\ can escape from these systems?",
            r"How does chaos affect\\ these escapes?"
        ).scale(1.8)
        for i in question:
            self.play(FadeIn(i, shift=UP))
            self.wait(0.2)
            self.next_section()
        self.play(FadeOut(question, shift=UP))
        self.wait(0.2)


class VariableScene(ThreeDScene):
    def construct(self):
        ax = ThreeDAxes()

        L = 15.
        r = 8.
        theta = 5 * PI / 6
        z = 3
        E = 15.

        newt = MNCenter()
        test = MNCTestParticle3D(L, E, r, z, theta, axis=ax)
        center = (test.get_center() - newt.get_center()) / 2
        newt.shift(-center)
        test.shift(-center)

        sphere_line = Line(
            newt.get_center(), test.get_center()
        )
        sphere_line.add_updater(
            lambda obj: obj.put_start_and_end_on(
                newt.get_center(), test.get_center()
            )
        )

        polar_line = Line(
            newt.get_center(), test.get_center() - (test.get_center() - newt.get_center())[2] * OUT
        )
        polar_line.add_updater(
            lambda obj: obj.put_start_and_end_on(
                newt.get_center(), test.get_center() - (test.get_center() - newt.get_center())[2] * OUT
            )
        )

        proj_line = Line(
            test.get_center(), test.get_center() - (test.get_center() - newt.get_center())[2] * OUT
        )
        proj_line.add_updater(
            lambda obj: obj.put_start_and_end_on(
                test.get_center(), test.get_center() - (test.get_center() - newt.get_center())[2] * OUT
            )
        )

        r_text = MathTex('r').scale(1.5).next_to(
            polar_line.get_center(), IN
        ).shift(UP)
        z_text = MathTex('z').scale(1.5).next_to(
            proj_line.get_center(), UP
        )
        R_text = MathTex(r'\sqrt{r^2+z^2}=R').scale(1.5).next_to(
            sphere_line.get_center(), (OUT + DOWN), buff=LARGE_BUFF
        )
        field = Arrow(
            *proj_line.get_start_and_end()[::-1]
        ).set_color(YELLOW).shift(2 * UP)
        B_text = MathTex('B').scale(1.5).set_color(YELLOW).next_to(field, UP)

        theta0 = Line(
            newt.get_center(), newt.get_center() + 2 * UP
        )
        angle = Angle(
            theta0.copy().shift(theta0.get_start()[2]*IN),
            polar_line.copy().shift(polar_line.get_start()[2]*IN),
            radius=1.5
        ).shift(theta0.get_start()[2]*OUT)
        theta_text = MathTex(r'\theta').scale(1.5).next_to(angle, UP).shift(0.125 * LEFT)

        text_1 = VGroup(r_text, z_text, theta_text, B_text, R_text)

        self.set_camera_orientation(**default_camera)
        self.add_fixed_orientation_mobjects(*text_1)
        self.remove(*text_1)
        self.play(GrowFromCenter(newt), GrowFromCenter(test))
        self.wait()
        self.play(
            *[Create(i) for i in (polar_line, proj_line, sphere_line, theta0, angle)],
            Write(text_1[:-2])
        )
        self.wait(0.2)
        self.next_section()
        self.play(
            Create(field), Write(text_1[-2])
        )
        self.wait(0.2)
        self.next_section()
        self.play(
            Write(text_1[-1])
        )
        self.wait(0.2)
        self.next_section()
        test_text = MathTex('(m, q)').scale(1.5).next_to(test, OUT, buff=MED_LARGE_BUFF).set_color(test.get_color())
        newt_text = MathTex('M').scale(1.5).next_to(newt, DOWN, buff=MED_LARGE_BUFF).set_color(newt.get_color())
        text_2 = VGroup(test_text, newt_text)
        self.add_fixed_orientation_mobjects(*text_2)
        self.remove(*text_2)
        self.play(
            Write(text_2[0])
        )
        self.wait(0.2)
        self.next_section()
        self.play(
            Write(text_2[1])
        )
        self.wait(0.2)
        self.next_section()
        self.play(
            *[Uncreate(i) for i in (polar_line, proj_line, field, sphere_line, angle, theta0, theta0, newt, test)],
            *[Unwrite(i) for i in (text_1, text_2)]
        )
        self.wait(0.2)


class HamiltonianScene(Scene):
    def construct(self):
        TitleScreen = Tex(r"Magnetized Newtonian\\ center").scale(2.5)
        self.play(Write(TitleScreen, run_time=0.5))
        self.wait(0.2)
        self.next_section()
        self.play(Unwrite(TitleScreen, run_time=0.5))

        lagrange_def = MathTex(r'\mathcal{L}=', 'T', '-', 'U')
        lagrange_ful = MathTex(r'\mathcal{L}=', r'\frac 12m(\dot{r}^2+r^2\dot{\theta}^2+\dot{z}^2)', '+',
                               r'\left(\frac12qBr^2\dot{\theta}+{GMm\over R}\right)',
                               )
        lagrange_ful[1][3].set_color(YELLOW)
        lagrange_ful[-1][4].set_color(YELLOW)
        lagrange_ful[-1][5].set_color(BLUE)
        lagrange_ful[-1][-4].set_color(YELLOW)
        lagrange_ful[-1][-5].set_color(RED)
        label1 = MathTex('(1)').next_to(lagrange_ful, RIGHT).to_edge(RIGHT)
        self.play(Write(lagrange_def, run_time=0.5))
        self.wait()
        self.play(ReplacementTransform(lagrange_def, lagrange_ful), Write(label1))
        self.wait(0.2)
        self.next_section()

        lagrange_less = MathTex(r'\mathcal{L}=\frac 12(\dot{r}^2+r^2\dot{\theta}^2+\dot{z}^2)+\frac 12 r^2\dot{\theta}+{1\over R}')
        label2 = MathTex('(2)').next_to(lagrange_less, RIGHT).to_edge(RIGHT)
        self.play(FadeOut(lagrange_ful, shift=UP), FadeOut(label1, shift=UP),
                  FadeIn(lagrange_less, shift=UP), FadeIn(label2, shift=UP))
        self.wait(0.2)
        self.next_section()

        char_l = MathTex(r"l_c=\sqrt[3]{GM t_c^2}", tex_to_color_map={"M": ORANGE}).to_corner(DL).shift(RIGHT * 1.8 + 2 * UP)
        label3 = MathTex('(3)').next_to(char_l, RIGHT, buff=1)
        char_t = MathTex(r"t_c=mq/B", tex_to_color_map={"mq": BLUE, "B": YELLOW}).next_to(label3, RIGHT, buff=2)
        label4 = MathTex('(4)').next_to(char_t, RIGHT, buff=1)
        self.play(*[FadeIn(i, shift=UP) for i in (char_l, char_t, label3, label4)])
        self.wait(0.2)
        self.next_section()

        ham = MathTex(
            r"h", r"=\mathcal{H}=\frac{1}{2}(p_r^2+p_z^2)+\frac{1}{2}\left({",
            "L", r"\over r}-\frac{r}{2}\right)-\frac{1}{R}"
        )
        ham.add(MathTex('(5)').next_to(ham, RIGHT).to_edge(RIGHT))

        pr = MathTex(r'\dot{p}_r={', r'L', r'\over r^3}-\frac{r}{4}-\frac{r}{R^3}')
        ptheta = MathTex(r'\dot{p}_\theta', '=0')
        pz = MathTex(r'\dot{p}_z=-\frac{z}{R^3}')

        p = VGroup(pr, ptheta, pz).arrange(DOWN, aligned_edge=LEFT).to_edge(LEFT, buff=1)
        label678 = VGroup(*[MathTex(i) for i in ('(6)', '(7)', '(8)')]).arrange(DOWN, buff=0.5).next_to(p, RIGHT,
                                                                                                        buff=1)
        leqn = MathTex(r'\Rightarrow p_\theta=L').next_to(p[1], RIGHT)
        r = MathTex(r'\dot{r}=p_r')
        theta = MathTex(r'\dot{\theta}={', r'L', r'\over r^2}-\frac 12')
        z = MathTex(r'\dot{z}=p_z')

        q = VGroup(r, theta, z).arrange(DOWN, aligned_edge=LEFT).to_edge(RIGHT, buff=3.25)
        label91011 = VGroup(*[MathTex(i) for i in ('(9)', '(10)', '(11)')]).arrange(DOWN, buff=0.5).next_to(q, RIGHT,
                                                                                                        buff=1)

        self.play(
            *[FadeOut(i, shift=UP) for i in (lagrange_less, label2)],
            FadeIn(ham, shift=UP)
        )
        self.wait(0.2)
        self.next_section()
        self.play(
            ham.animate.to_edge(UP),
            *[i.animate.shift(1.5*DOWN) for i in (char_l, char_t, label3, label4)]
        )
        self.play(*[Write(i) for i in (p, q, label678, label91011)], run_time=1)
        self.wait(0.2)
        self.next_section()
        ham[0].set_color(GREEN)
        self.play(Indicate(ham[0], color=GREEN, scale_factor=1.5), run_time=1.5)
        self.wait(0.2)
        self.next_section()
        ham[2].set_color(GREEN)
        self.play(Indicate(ham[2], color=GREEN, scale_factor=1.5), run_time=1.5)
        self.wait(0.2)
        self.next_section()
        self.play(Indicate(p[1], color=GREEN, scale_factor=1.5), run_time=1.5)
        self.wait(0.2)
        self.next_section()
        self.play(Write(leqn))
        self.wait(0.2)
        self.next_section()
        print(self.mobjects)
        all_text = VGroup()
        for i in self.mobjects:
            print(i)
            if isinstance(i, (MathTex, VGroup)):
                all_text.add(i)
        self.play(FadeOut(all_text))
        self.wait(0.2)


class NoTheta(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(**default_camera)
        L = -2.
        r = 2.
        z = 0.5
        E = 0.3
        theta = PI/2
        axis = ThreeDAxes()
        newt = MNCenter()
        test = MNCTestParticle3D(L, E, r, z, theta=theta, axis=axis, speed=2)

        def camera_move(obj):
            self.set_camera_orientation(theta=(test.qp[0, 1]-90*DEGREES))

        self.play(Create(axis))
        self.play(GrowFromCenter(test), GrowFromCenter(newt))
        self.move_camera(phi=90 * DEGREES, theta=(test.qp[0, 1] - 90 * DEGREES))
        self.add_updater(camera_move)
        test.move()

        def get_trajectory():
            pos = cartesian_to_cylindrical(test.get_center())
            return np.array([pos[0], pos[2], 0], float)

        path = TracedPath(
            get_trajectory, stroke_color=test.get_color(), stroke_width=DEFAULT_STROKE_WIDTH, stroke_opacity=1.0
        )
        self.add_fixed_in_frame_mobjects(path)
        self.wait(8)
        path.clear_updaters()
        test.stop()
        self.play(Uncreate(test), Uncreate(newt), Uncreate(axis))
        ax = Axes(x_length=(config.frame_height + 2.5)/2, x_range=(0, 7, 1))
        self.add_fixed_in_frame_mobjects(ax)
        ax.shift(-ax.c2p(*ORIGIN))
        ax_label = ax.get_axis_labels(x_label=MathTex('r'), y_label=MathTex('z'))
        self.add_fixed_in_frame_mobjects(ax_label)
        self.remove(ax, ax_label)
        self.play(Create(ax), Write(ax_label))
        self.wait(0.2)
        self.next_section()
        self.play(Uncreate(ax), Uncreate(path), Unwrite(ax_label))
        self.wait(0.2)


class ZVCScene(Scene):
    def construct(self):
        ham = MathTex(
            r"h", "=", r"\mathcal{H}=\frac{1}{2}(p_r^2+p_z^2)+", r"\frac{1}{2}\left({",
            "L", r"\over r}-\frac{r}{2}\right)-\frac{1}{R}"
        )
        label1 = MathTex('(1)').next_to(ham, RIGHT, buff=1)
        ZVC_eqn = MathTex(
            r"h", "=", r"\frac{1}{2}\left({",
            "L", r"\over r}-\frac{r}{2}\right)-\frac{1}{R}"
        )
        label12 = MathTex('(12)').move_to(label1.get_center())
        self.play(Write(ham, run_time=0.5), Write(label1, run_time=0.5))
        self.wait(0.2)
        self.next_section()
        ZVC_text = Tex("Zero Velocity Curve (ZVC)").to_edge(UR)
        self.play(
            FadeOut(ham[2]),
            *[i.animate.move_to(j) for i, j in [
                [ham[0], ZVC_eqn[0]],
                [ham[1], ZVC_eqn[1]],
                [ham[3], ZVC_eqn[2]],
                [ham[4], ZVC_eqn[3]],
                [ham[5], ZVC_eqn[4]]
            ]],
            ReplacementTransform(label1, label12),
            Write(ZVC_text, run_time=0.5)
        )
        self.wait(0.2)
        self.next_section()
        stretch = 2
        ax = Axes(
            x_length=3 * stretch,
            x_range=(0, 3)
        ).shift(RIGHT)
        ax_label = ax.get_axis_labels(x_label=MathTex('r'), y_label=MathTex('z'))
        E = ValueTracker(0.6)
        L = ValueTracker(-1.)
        zvc = MNCZeroVelocityCurve(L, E, ax)

        text_h = DecimalNumber(
            energy_E_to_h(L.get_value(), E.get_value())[0],
            num_decimal_places=3
        ).move_to(4 * LEFT + 0.25 * UP)
        text_L = DecimalNumber(
            L.get_value(),
            num_decimal_places=3
        ).next_to(text_h, DOWN)
        text = VGroup(MathTex('h', '=').next_to(text_h, LEFT), text_h)
        text1 = VGroup(MathTex('L', '=').next_to(text[0], DOWN), text_L)
        text = VGroup(text, text1)

        text_h.add_updater(lambda d: d.set_value(energy_E_to_h(L.get_value(), E.get_value())[0]))
        text_L.add_updater(lambda d: d.set_value(L.get_value()))

        self.play(
            *[i.animate.move_to(j) for i, j in [
                [ham[0], text[0][0][0]],
                [ham[4], text[1][0][0]]
            ]],
            *[FadeOut(i) for i in (ham[1], ham[3], ham[5], label12)],
            *[Write(i) for i in (ax_label, text[0][0][1], text[0][1], text[1][0][1], text[1][1])],
            *[Create(i) for i in (ax, zvc)]
        )
        self.remove(ZVC_eqn)
        self.add(ax, ax_label, zvc, text)
        self.wait(0.2)
        self.next_section()
        particle = MNCTestParticle2D(L.get_value(), E.get_value(), 1.0, 0, axis=ax)
        particle.speed.set_value(5.0)
        particle2 = MNCTestParticle2D(L.get_value(), E.get_value(), 1.8, 0.2, axis=ax)
        particle2.speed.set_value(5.0)
        particle3 = MNCTestParticle2D(L.get_value(), E.get_value(), 1.5, -2, axis=ax)
        particle3.speed.set_value(5.0)
        particle.move()
        particle2.move()
        particle3.move()
        self.add(particle, particle2, particle3)
        self.wait(10)
        self.play(Uncreate(particle), Uncreate(particle2), Uncreate(particle3))
        self.next_section()
        self.play(E.animate.set_value(0.4))
        self.play(L.animate.set_value(-0.5))
        self.play(E.animate.set_value(0.6))
        self.play(L.animate.set_value(-1))
        self.wait(0.2)
        self.next_section()

        def vector_function(pos):
            r, z = ax.p2c(pos)
            return RIGHT*(L.get_value()**2/r**3-r/4-r/(r**2+z**2)**1.5)-UP*z/(r**2+z**2)**1.5
        vector_field = ArrowVectorField(
            vector_function,
            x_range=[ax.c2p(*(0.2*RIGHT))[0], ax.c2p(*(2.5*RIGHT))[0], 0.4],
            y_range=[-ax.c2p(*(3*UP))[1], ax.c2p(*(3*UP))[1], 0.4],
            color=WHITE
        )
        self.play(Create(vector_field))
        self.wait(0.2)
        self.next_section()

        r_eq = MNCValues(L.get_value(), E.get_value()).get_r_eq()
        minima_point = Dot(ax.c2p(r_eq, 0))
        self.play(E.animate.set_value(0.005), run_time=2)
        zvc.stop()
        self.play(ReplacementTransform(zvc, minima_point), E.animate.set_value(0.), Uncreate(vector_field))
        self.wait(0.2)


class CircularOrbitScene(ThreeDScene):
    def construct(self):
        stretch = 2
        ax = Axes(
            x_length=3 * stretch,
            x_range=(0, 3)
        ).shift(RIGHT)
        ax_label = ax.get_axis_labels(x_label=MathTex('r'), y_label=MathTex('z'))
        E_old = ValueTracker(0.)
        L = ValueTracker(-1.)

        text_h = DecimalNumber(
            energy_E_to_h(L.get_value(), E_old.get_value())[0],
            num_decimal_places=3
        ).move_to(4 * LEFT + 0.25 * UP)
        text_L = DecimalNumber(
            L.get_value(),
            num_decimal_places=3
        ).next_to(text_h, DOWN)
        text1 = VGroup(MathTex('h', '=').next_to(text_h, LEFT), text_h)
        text2 = VGroup(MathTex('L', '=').next_to(text1[0], DOWN), text_L)
        text = VGroup(text1, text2)

        text_h.add_updater(lambda d: d.set_value(energy_E_to_h(L.get_value(), E_old.get_value())[0]))
        text_L.add_updater(lambda d: d.set_value(L.get_value()))

        r_eq = MNCValues(L.get_value(), E_old.get_value()).get_r_eq()
        minima_point = Dot(ax.c2p(r_eq, 0))
        ZVC_text = Tex("Zero Velocity Curve (ZVC)").to_edge(UR)

        self.add(ax_label, text, ax, minima_point, ZVC_text)
        phi, _, _, _, _ = self.camera.get_value_trackers()
        self.play(FadeOut(text, ZVC_text), minima_point.animate.set_color(BLUE), phi.animate.set_value(-30*DEGREES), run_time=0.5)
        path = TracedPath(
            minima_point.get_center,
            stroke_color=minima_point.get_color(),
            stroke_width=DEFAULT_STROKE_WIDTH,
            stroke_opacity=1.0
        )
        self.add(path)
        self.play(
            Rotate(minima_point, axis=DOWN, about_point=ax.c2p(*ORIGIN), rate_func=linear, run_time=7, angle=6*PI)
        )
        self.play(
            FadeIn(text), Uncreate(path), minima_point.animate.set_color(WHITE), phi.animate.set_value(0 * DEGREES),
            run_time=0.5
        )
        self.wait(0.2)
        self.next_section()
        coord = MathTex(r'(r_{\text{eq}}, 0)').next_to(minima_point, UP)
        quart = MathTex(r'(r_{\text{eq}})^4+4(r_{\text{eq}})-4L^2=0', ).next_to(minima_point, DOWN).shift(RIGHT*1.5)
        label13 = MathTex('(13)').next_to(quart, RIGHT, buff=1)
        self.play(Write(coord))
        self.play(*[Write(i) for i in (quart, label13)], run_time=0.5)
        self.wait(0.2)
        self.next_section()
        U_eq = MathTex("U_{eq}").next_to(text1[0][1], RIGHT)
        text_h.clear_updaters()
        self.remove(text_h)
        self.play(Indicate(U_eq, color=GREEN, scale_factor=1.5), run_time=1.5)
        self.wait(0.2)
        self.next_section()
        E = ValueTracker(0.005)
        zvc = MNCZeroVelocityCurve(L, E, ax)
        text_h.add_updater(lambda d: d.set_value(energy_E_to_h(L.get_value(), E.get_value())[0]))
        self.play(
            *[Unwrite(i, run_time=0.5) for i in (coord, quart, label13)],
            ReplacementTransform(minima_point, zvc),
            FadeOut(U_eq), FadeIn(text_h)
        )
        self.play(E.animate.set_value(0.99))
        self.wait(0.2)
        self.next_section()
        asymptote = MNCAsymptotes(L, E, ax)
        self.add(asymptote)
        self.play(E.animate.set_value(1.1), run_time=3)
        self.wait(0.2)
        self.next_section()
        self.play(E.animate.set_value(1.05), run_time=0.5)
        self.play(E.animate.set_value(1.15), run_time=0.5)
        self.play(E.animate.set_value(1.1), run_time=0.5)
        self.wait(0.2)
        self.next_section()
        self.play(E.animate.set_value(1.2))
        self.wait(0.2)
        self.next_section()
        self.play(E.animate.set_value(1.0))
        U_esc = MathTex("U_{esc, min}").next_to(text1[0][1], RIGHT)
        text_h.clear_updaters()
        self.remove(text_h)
        self.play(Indicate(U_esc, color=GREEN, scale_factor=1.5), run_time=1.5)
        self.wait(0.2)
        self.next_section()
        self.play(FadeOut(U_esc), FadeIn(text_h))
        text_h.add_updater(lambda d: d.set_value(energy_E_to_h(L.get_value(), E_old.get_value())[0]))
        self.play(E.animate.set_value(1.2))
        particle = MNCTestParticle2D(L.get_value(), E.get_value(), 1.0, 0, axis=ax)
        particle.speed.set_value(1.5)
        particle.move()
        path = TracedPath(
            particle.get_center, stroke_color=particle.get_color(), stroke_width=DEFAULT_STROKE_WIDTH, stroke_opacity=1.0
        )
        self.add(particle, path)
        self.wait(5)


class RescaleScene(Scene):
    def construct(self):
        Umin_text = Tex("Equilibrium energy").scale(1.5)
        Umin_eqn = MathTex(
            "U_{eq}", "=\\frac{1}{2}\\left(\\frac{L}{r_{eq}}-\\frac{r_{eq}}{2}\\right)-\\frac{1}{r_{eq}}"
        )
        rmin_eqn = MathTex("r_{eq}^4+4r_{eq}-4L^2=0")
        Umin = VGroup(Umin_text, Umin_eqn, rmin_eqn).arrange(DOWN).to_edge(UP)
        Umin[2].shift(0.5*LEFT)
        Uesc_eqn = MathTex(
            "U_{esc, min}", "=\\frac{1}{2}(|L|-L)"
        )
        Uesc_text = Tex("Energy of escape").scale(1.5)
        Uesc = VGroup(Uesc_text, Uesc_eqn).arrange(DOWN).next_to(Umin, DOWN, buff=LARGE_BUFF)
        label14 = MathTex('(14)').to_edge(RIGHT)
        Umin_eqn.add(label14.move_to(label14.get_x()*RIGHT+Umin_eqn.get_y()*UP))
        label15 = MathTex('(15)')
        rmin_eqn.add(label15.move_to(label14.get_x()*RIGHT+rmin_eqn.get_y()*UP))
        label16 = MathTex('(16)')
        Uesc_eqn.add(label16.move_to(label14.get_x()*RIGHT+Uesc_eqn.get_y()*UP))

        self.play(Write(Umin, run_time=0.5))
        self.wait(0.2)
        self.next_section()
        self.play(Write(Uesc, run_time=0.5))
        self.wait(0.2)
        self.next_section()
        self.play(*[Unwrite(i, run_time=0.5) for i in (Uesc, Umin)])
        text_1 = Tex("Define a new energy scale ", "E").scale(1.2).to_edge(UP)
        E0 = MathTex("h=U_{eq} \\Rightarrow E=0")
        E1 = MathTex("h=U_{esc,min} \\Rightarrow E=1")
        rescale = MathTex("h=", r"U_{eq}", "+", "E", "(", r"U_{esc,min}", "-", r"U_{eq}", ")").scale(1.1)
        text_2 = Tex("*For a given (E, $r_0$, $z_0$, $p_{r0}$), trajectories are").scale(0.9)
        text_3 = Tex("symmetrical under a parity transformation of L").scale(0.9)
        rescale_group = VGroup(
            text_1, E0, E1, rescale, text_2, text_3
        ).arrange(DOWN).scale(1.5).to_edge(UP)
        label19 = MathTex("(19)").scale(1.5).to_edge(RIGHT)
        rescale.add(label19.move_to(label19.get_x()*RIGHT+rescale.get_y()*UP))
        E0.add(MathTex("(17)").scale(1.5).move_to(label19.get_x()*RIGHT+E0.get_y()*UP))
        E1.add(MathTex("(18)").scale(1.5).move_to(label19.get_x()*RIGHT+E1.get_y()*UP))
        text_2.shift(0.5*DOWN)
        text_3.shift(0.5*DOWN)
        for i in rescale_group:
            self.play(Write(i, run_time=0.5))
            self.wait(0.2)
            self.next_section()
        self.play(Unwrite(rescale_group, run_time=0.5))
        self.wait(0.2)


class SymmetryScene(Scene):
    def construct(self):
        stretch = 1
        ax = Axes(
            x_length=3 * stretch,
            x_range=(0, 3)
        ).shift(RIGHT)
        ax_label = ax.get_axis_labels(x_label=MathTex('r'), y_label=MathTex('z'))
        L = ValueTracker(1.)
        ax_neg = Axes(
            x_length=3 * stretch,
            x_range=(0, 3)
        ).shift(RIGHT)
        ax_label_neg = ax.get_axis_labels(x_label=MathTex('r'), y_label=MathTex('z'))
        L_neg = ValueTracker(-1.)
        E = ValueTracker(0.6)
        zvc = MNCZeroVelocityCurve(L, E, ax)
        zvc_neg = MNCZeroVelocityCurve(L_neg, E, ax_neg)

        text_E = DecimalNumber(
            E.get_value(),
            num_decimal_places=3
        ).next_to(ax, LEFT).shift(0.25 * UP)
        text_L = DecimalNumber(
            L.get_value(),
            num_decimal_places=3
        ).next_to(text_E, DOWN)
        text = VGroup(MathTex('E', '=').next_to(text_E, LEFT), text_E)
        text1 = VGroup(MathTex('L', '=').next_to(text[0], DOWN), text_L)
        text = VGroup(text, text1)

        text_E.add_updater(lambda d: d.set_value(E.get_value()))
        text_L.add_updater(lambda d: d.set_value(L.get_value()))

        text_E_neg = DecimalNumber(
            E.get_value(),
            num_decimal_places=3
        ).next_to(ax_neg, LEFT).shift(0.25 * UP)
        text_L_neg = DecimalNumber(
            L_neg.get_value(),
            num_decimal_places=3
        ).next_to(text_E, DOWN)
        text_neg = VGroup(MathTex('E', '=').next_to(text_E_neg, LEFT), text_E_neg)
        text1_neg = VGroup(MathTex('L', '=').next_to(text_neg[0], DOWN), text_L_neg)
        text_neg = VGroup(text_neg, text1_neg)

        text_E_neg.add_updater(lambda d: d.set_value(E.get_value()))
        text_L_neg.add_updater(lambda d: d.set_value(L_neg.get_value()))
        group = VGroup(ax, ax_label, text, zvc).to_edge(LEFT)
        group_neg = VGroup(ax_neg, ax_label_neg, text_neg, zvc_neg).to_edge(RIGHT)
        r = ValueTracker(1.0)
        T = ValueTracker(10)
        y_max = config["frame_y_radius"]
        divider = Line(start=UP*y_max, end=DOWN*y_max)
        self.play(*[FadeIn(i) for i in (divider, group[:-1], group_neg[:-1])])
        self.play(Create(group[-1]), Create(group_neg[-1]))
        self.wait(0.2)
        self.next_section()
        r_eq = MNCValues(L.get_value(), E.get_value()).get_r_eq()
        self.play(E.animate.set_value(0.01), run_time=1.0)
        asymptotes = MNCAsymptotes(L, E, ax)
        asymptotes_neg = MNCAsymptotes(L_neg, E, ax_neg)
        self.add(asymptotes, asymptotes_neg)
        self.play(E.animate.set_value(1.0), run_time=1.0)
        self.wait(0.2)
        self.next_section()
        self.play(E.animate.set_value(1.1))
        self.wait(0.2)
        self.next_section()
        self.play(E.animate.set_value(0.6))
        self.wait(0.2)
        self.next_section()
        trajectory = MNCTrajectory2D(L, E, r, ax, T)
        trajectory_neg = MNCTrajectory2D(L_neg, E, r, ax_neg, T)
        self.play(*[Create(i, rate_func=linear) for i in (trajectory, trajectory_neg)])
        self.wait(0.2)
        self.next_section()
        self.play(E.animate.set_value(0.4), run_time=2.0)
        self.play(L.animate.set_value(0.5), L_neg.animate.set_value(-0.5), r.animate.set_value(0.5), run_time=2.0)
        self.play(E.animate.set_value(1.1), run_time=2.0)
        self.play(L.animate.set_value(1), L_neg.animate.set_value(-1), r.animate.set_value(1), run_time=2.0)
        self.play(E.animate.set_value(0.6), run_time=2.0)
        self.wait(0.2)


class EscapeDefinitionScene(Scene):
    def construct(self):
        slide_title = Tex("Escape Dynamics").scale(3)
        what = VGroup(
            Tex("How do we qualify a"),
            Tex("trajectory as escaping?")
        ).arrange(DOWN).scale(2)
        title_group = VGroup(slide_title, what).arrange(DOWN)
        self.play(Write(title_group[0], run_time=0.5))
        self.wait(0.2)
        self.next_section()
        self.play(Write(title_group[1], run_time=0.5))
        self.wait(0.2)
        self.next_section()
        self.play(Unwrite(title_group, run_time=0.5))

        stretch = 3
        ax = Axes(
            x_length=3 * stretch,
            x_range=(0, 3),
            y_range=(0, 13)
        )
        ax_label = ax.get_axis_labels(x_label=MathTex('r'), y_label=MathTex('z'))
        L = ValueTracker(1.0)
        E = ValueTracker(1.05)
        r_esc = 1.0
        pr_esc = 0.2
        r_inesc = 0.8475439238215814
        pr_inesc = -0.040705422983666084
        r = ValueTracker(r_esc)
        pr = ValueTracker(pr_esc)
        T = ValueTracker(70)
        trajectory = MNCTrajectory2D(L, E, r, ax, T, pr)
        trajectory_plot = VGroup(ax, ax_label, trajectory)
        self.play(Create(trajectory_plot))
        self.wait(0.2)
        self.next_section()
        new_ax = Axes(
            x_length=3 * stretch,
            x_range=(0, 3),
            y_range=(0, 20)
        )
        new_trajectory = MNCTrajectory2D(L, E, r, new_ax, T, pr)
        self.play(
            ReplacementTransform(ax, new_ax),
            ReplacementTransform(trajectory, new_trajectory)
        )
        self.wait(0.2)
        self.next_section()
        self.play(
            r.animate.set_value(r_inesc),
            pr.animate.set_value(pr_inesc)
        )
        self.wait(0.2)
        self.next_section()

        ZVC = MNCZeroVelocityCurve(L, E, new_ax)
        amplitudes = Amplitudes(new_trajectory, new_ax)
        asymptotes = MNCAsymptotes(L, E, new_ax)
        is_escaping = Text("escaping").to_corner(UR)

        def definition_updater(mob):
            if amplitudes.get_right()[0] > asymptotes.get_right()[0]:
                mob.become(Text("inescaping").to_corner(UR))
            elif amplitudes.get_right()[0] < asymptotes.get_right()[0]:
                mob.become(Text("escaping").to_corner(UR))

        is_escaping.add_updater(definition_updater)
        self.play(Create(ZVC))
        self.wait(0.2)
        self.next_section()

        def vector_function(pos):
            r, z = ax.p2c(pos)
            return RIGHT*(L.get_value()**2/r**3-r/4-r/(r**2+z**2)**1.5)-UP*z/(r**2+z**2)**1.5
        vector_field = ArrowVectorField(
            vector_function,
            x_range=[ax.c2p(*(0.2*RIGHT))[0], ax.c2p(*(3*RIGHT))[0], 0.4],
            y_range=[ax.c2p(*(0*UP))[1], ax.c2p(*(13*UP))[1], 0.4],
            color=WHITE
        )
        self.play(FadeIn(vector_field))
        self.wait(0.2)
        self.next_section()
        self.play(Create(asymptotes), FadeOut(vector_field))
        self.wait(0.2)
        self.next_section()
        self.play(Create(amplitudes))
        self.wait(0.2)
        self.next_section()
        self.play(Write(is_escaping))
        self.wait(0.2)
        self.next_section()
        self.play(
            r.animate.set_value(r_esc),
            pr.animate.set_value(pr_esc)
        )
        self.wait(0.2)
        self.next_section()
        self.play(
            r.animate.set_value(r_inesc),
            pr.animate.set_value(pr_inesc)
        )
        self.wait(0.2)
        self.next_section()
        self.play(T.animate.set_value(160))
        self.wait(0.2)
        self.next_section()
        self.play(*[FadeOut(i) for i in (ZVC, asymptotes, amplitudes, is_escaping, trajectory_plot)])
        sentence = VGroup(
            Tex("An MNC trajectory escapes when its"),
            Tex("amplitudes of oscillation"),
            Tex("are inbetween the"),
            Tex("asymptotes of escape.")
        ).arrange(DOWN).scale(1.5)
        sentence[1].set_color(GREEN)
        sentence[3].set_color(RED)
        self.play(Write(sentence, run_time=0.5))
        self.wait(0.2)
        self.next_section()
        self.play(Unwrite(sentence, run_time=0.5))
        self.wait(0.2)


class EscapeQuantities(Scene):
    def construct(self):
        escape_title = Tex("Escape Quantities").scale(3)
        quantity_list = BulletedList(
            r"Escape Time",
            r"Time it takes for a\\trajectory to escape",
            r"Escape Pass",
            r"How many times a\\trajectory crosses the\\equatorial plane $z=0$",
        )
        self.play(Write(escape_title, run_time=0.5))
        self.wait(0.5)
        self.play(escape_title.animate.scale(1/2).to_corner(UL))
        quantity_list.next_to(escape_title, DOWN)
        for i in (1, 3):
            quantity_list[i].shift(0.2*RIGHT)
        stretch = 2
        ax = Axes(
            x_length=3 * stretch,
            x_range=(0, 3),
            y_range=(-4, 14)
        ).to_corner(UR)
        ax_label = ax.get_axis_labels(x_label=MathTex('r'), y_label=MathTex('z'))
        L = 1.0
        E = 1.05
        r = 2.5485392707361787
        pr = -0.12211626895099847
        T = 43.79
        speed = ValueTracker(T/10)
        zvc = MNCZeroVelocityCurve(ValueTracker(L), ValueTracker(E), ax)
        asymptotes = MNCAsymptotes(ValueTracker(L), ValueTracker(E), ax)
        plot = VGroup(ax, ax_label, zvc, asymptotes)

        particle = MNCTestParticle2D(L, E, r, pr=pr, speed=speed, axis=ax)
        text_time = DecimalNumber(
            0.,
            num_decimal_places=2
        )
        text = VGroup(Tex('Escape Time', ' ='), text_time).arrange(RIGHT).next_to(ax, DOWN)

        self.play(
            Write(quantity_list[:2], run_time=0.5),
            Create(plot), Write(text, run_time=0.5), Create(particle)
        )
        self.wait(0.2)
        self.next_section()
        time_trajectory = TracedPath(
            particle.get_center, stroke_color=particle.get_color(), stroke_width=DEFAULT_STROKE_WIDTH
        )

        def text_time_updater(mob, dt):
            mob.set_value(mob.get_value()+speed.get_value()*dt)
        self.add(time_trajectory)
        particle.move()
        text_time.add_updater(text_time_updater)
        self.wait(10.0)
        particle.stop()
        text_time.clear_updaters()
        time_trajectory.clear_updaters()
        self.wait(0.2)
        self.next_section()

        pass_particle = MNCTestParticle2D(L, E, r, pr=pr, speed=speed, axis=ax)
        text_pass = DecimalNumber(
            0.,
            num_decimal_places=0
        )
        text_pass_full = VGroup(Tex('Escape Pass', ' ='), text_pass).arrange(RIGHT).next_to(ax, DOWN)
        pass_trajectory = TracedPath(
            pass_particle.get_center, stroke_color=pass_particle.get_color(), stroke_width=DEFAULT_STROKE_WIDTH
        )

        def text_pass_updater(mob):
            z_sign = ax.p2c(pass_particle.points[-1])[1]*ax.p2c(pass_particle.points[-3])[1]
            print(ax.p2c(pass_particle.points[-1]), z_sign)
            if z_sign < 0:
                mob.set_value(mob.get_value()+1)

        self.play(
            Unwrite(text, run_time=0.5),
            Uncreate(particle), Uncreate(time_trajectory)
        )
        self.play(Write(quantity_list[2:], run_time=0.5), Write(text_pass_full, run_time=0.5), Create(pass_particle))
        self.add(pass_trajectory)
        self.wait(0.2)
        self.next_section()
        pass_particle.move()
        self.wait(0.2)
        text_pass.add_updater(text_pass_updater)
        self.wait(9.8)
        pass_particle.stop()
        text_pass.clear_updaters()
        self.wait(0.2)


class EscapePlot(Scene):
    def construct(self):
        escape_time = ImageMobject("escape_time.png").scale(0.8)
        escape_time_scale = ImageMobject("escape_time_scale.png").scale(0.8)
        escape_time_text = Tex(r"Fig 1: Escape time plot\\ (L, E)=(1, 1.05)").next_to(escape_time, DOWN)
        escape_axis = Axes(
            x_range=(0.4187087820553533, 2.6066446574056354),
            y_range=(-1.361310153581102, 1.361310153581102),
            x_length=escape_time_scale.width,
            y_length=escape_time_scale.height
        ).shift(
            0.06 * UP * escape_time_scale.height +
            0.03 * LEFT * escape_time_scale.width
        )
        escape_time_group = Group(escape_axis, escape_time, escape_time_text).to_edge(LEFT)
        self.play(FadeIn(escape_time_group[1:]))
        self.wait(0.2)
        self.next_section()
        stretch = 2
        traj_axis = Axes(
            x_length=3 * stretch,
            x_range=(0, 3),
            y_range=(-5, 5)
        )
        ax_label = traj_axis.get_axis_labels(x_label=MathTex('r'), y_label=MathTex('z'))
        L = ValueTracker(1.0)
        E = ValueTracker(1.05)
        r = ValueTracker(2.0)
        pr = ValueTracker(0)
        T = ValueTracker(100)
        zvc = MNCZeroVelocityCurve(L, E, traj_axis)
        trajectory = MNCTrajectory2D(L, E, r, traj_axis, T, pr)
        asymptotes = MNCAsymptotes(L, E, traj_axis)
        trajectory_plot = VGroup(traj_axis, ax_label, trajectory, zvc, asymptotes).to_edge(RIGHT)
        ic_point = VGroup(
            Dot(escape_axis.c2p(r.get_value(), pr.get_value()), stroke_width=3, color=BLACK, stroke_color=WHITE),
            Circle(radius=0.5, stroke_opacity=0.5, color=BLACK).move_to(escape_axis.c2p(r.get_value(), pr.get_value()))
        )
        ic_point.add_updater(
            lambda mob: mob.move_to(escape_axis.c2p(r.get_value(), pr.get_value()))
        )
        self.play(Create(trajectory_plot), Create(ic_point))
        self.wait(0.2)
        self.next_section()
        self.play(r.animate.set_value(1.0), pr.animate.set_value(-0.3))
        self.wait(0.2)
        self.next_section()
        self.play(r.animate.set_value(1.05))
        self.play(r.animate.set_value(0.95))
        self.play(r.animate.set_value(1.1))
        self.wait(0.2)
        self.next_section()
        amplitudes = Amplitudes(trajectory, traj_axis)
        self.play(
            r.animate.set_value(1.1), pr.animate.set_value(0.1)
        )
        self.play(Create(amplitudes))
        trajectory_plot.add(amplitudes)
        self.wait(0.2)
        self.next_section()
        self.play(r.animate.set_value(1.15))
        self.play(r.animate.set_value(1.05))
        self.play(r.animate.set_value(1.15))
        self.wait(0.2)
        self.next_section()
        self.play(r.animate.set_value(1.5), pr.animate.set_value(0.3))
        self.wait(0.2)
        self.next_section()
        escape_time_zoomed = ImageMobject("escape_time_zoomed.png").scale(0.8)
        escape_time_zoomed_scale = ImageMobject("escape_time_zoomed_scale.png").scale(0.8)
        escape_time_zoomed_text = Tex(r"Fig 2: Zoomed escape time\\ plot (L, E)=(1, 1.05)").next_to(escape_time_zoomed, DOWN)
        escape_axis_zoomed = Axes(
            x_range=(1.4, 1.6),
            y_range=(0.2, 0.4),
            x_length=escape_time_zoomed_scale.width,
            y_length=escape_time_zoomed_scale.height
        ).shift(
            0.06 * UP * escape_time_zoomed_scale.height +
            0.035 * LEFT * escape_time_zoomed_scale.width
        )
        escape_time_zoomed_group = Group(escape_axis_zoomed, escape_time_zoomed, escape_time_zoomed_text).to_edge(LEFT)
        ic_point.clear_updaters().set_z_index(escape_time_zoomed_group.z_index+1)
        self.play(
            FadeOut(escape_time_group),
            FadeIn(escape_time_zoomed_group),
            ic_point.animate.move_to(escape_axis_zoomed.c2p(r.get_value(), pr.get_value()))
        )
        ic_point.add_updater(
            lambda mob: mob.move_to(escape_axis_zoomed.c2p(r.get_value(), pr.get_value()))
        )
        self.wait(0.2)
        self.next_section()
        self.play(r.animate.set_value(1.4))
        self.wait(0.2)
        self.next_section()
        self.play(r.animate.set_value(1.55))
        self.wait(0.2)


class EscapePlot1(Scene):
    def construct(self):
        L = ValueTracker(1.0)
        E = ValueTracker(1.05)
        r = ValueTracker(1.55)
        pr = ValueTracker(0.3)
        T = ValueTracker(250)

        escape_time_zoomed = ImageMobject("escape_time_zoomed.png").scale(0.8)
        escape_time_zoomed_scale = ImageMobject("escape_time_zoomed_scale.png").scale(0.8)
        escape_time_zoomed_text = Tex(r"Fig 2: Zoomed escape time\\ plot (L, E)=(1, 1.05)").next_to(escape_time_zoomed, DOWN)
        escape_axis_zoomed = Axes(
            x_range=(1.4, 1.6),
            y_range=(0.2, 0.4),
            x_length=escape_time_zoomed_scale.width,
            y_length=escape_time_zoomed_scale.height
        ).shift(
            0.06 * UP * escape_time_zoomed_scale.height +
            0.035 * LEFT * escape_time_zoomed_scale.width
        )
        escape_time_zoomed_group = Group(escape_axis_zoomed, escape_time_zoomed, escape_time_zoomed_text).to_edge(LEFT)
        ic_point = VGroup(
            Dot(escape_axis_zoomed.c2p(r.get_value(), pr.get_value()), stroke_width=3, color=BLACK, stroke_color=WHITE),
            Circle(radius=0.5, stroke_opacity=0.5, color=BLACK).move_to(escape_axis_zoomed.c2p(r.get_value(), pr.get_value()))
        )
        ic_point.add_updater(
            lambda mob: mob.move_to(escape_axis_zoomed.c2p(r.get_value(), pr.get_value()))
        )

        stretch = 2
        traj_axis_zoomed = Axes(
            x_length=3 * stretch,
            x_range=(0, 3),
            y_range=(0, 60, 10)
        )
        ax_label_zoomed = traj_axis_zoomed.get_axis_labels(x_label=MathTex('r'), y_label=MathTex('z'))
        zvc_zoomed = MNCZeroVelocityCurve(L, E, traj_axis_zoomed)
        trajectory_zoomed = MNCTrajectory2D(L, E, r, traj_axis_zoomed, T, pr)
        asymptotes_zoomed = MNCAsymptotes(L, E, traj_axis_zoomed)
        amplitudes_zoomed = Amplitudes(trajectory_zoomed, traj_axis_zoomed)
        trajectory_plot_zoomed = VGroup(
            traj_axis_zoomed, ax_label_zoomed, trajectory_zoomed, zvc_zoomed, asymptotes_zoomed, amplitudes_zoomed
        ).to_edge(RIGHT)
        self.add(trajectory_plot_zoomed, escape_time_zoomed_group, ic_point)
        self.wait(0.2)
        self.next_section()
        self.play(r.animate.set_value(1.5), T.animate.set_value(600), run_time=3)
        self.wait(0.2)


class EscapePlot2(Scene):
    def construct(self):
        escape_time = ImageMobject("escape_time.png").scale(0.8)
        escape_time_text = Tex(r"Fig 1: Escape time plot\\ (L, E)=(1, 1.05)").next_to(escape_time, DOWN)
        escape_time_group = Group(escape_time, escape_time_text).to_edge(LEFT)

        escape_pass = ImageMobject("escape_pass.png").scale(0.8)
        escape_pass_text = Tex(r"Fig 3: Escape pass plot\\ (L, E)=(1, 1.05)").next_to(escape_pass, DOWN)
        escape_pass_group = Group(escape_pass, escape_pass_text).to_edge(RIGHT)

        list = BulletedList(
            r"For E greater than 1, trajectories\\can escape along $z\pm\infty$",
            r"There exists initial conditions\\with escapes that are\\insensitive to perturbations",
            r"There exists initial conditions\\with escapes that are\\sensitive to perturbations",
            r"Initial conditions on boundaries\\of the two aformentioned regions\\have high escape times"
        ).scale(0.9).to_edge(RIGHT)
        self.add(escape_time_group)
        self.play(FadeIn(escape_pass_group))
        self.wait(0.2)
        self.next_section()
        self.play(FadeOut(escape_pass_group))
        for i in list:
            self.play(Write(i, run_time=0.5))
            self.wait(0.2)
            self.next_section()
        self.play(Unwrite(list, run_time=0.5), FadeOut(escape_time_group))
        self.wait(0.2)


class TrapDemo(Scene):
    def construct(self):
        trap_1 = ImageMobject("trap_demo_1.png").scale(0.5)
        trap_1.add(Tex("Fig 4a: L=2").next_to(trap_1, DOWN, buff=SMALL_BUFF).scale(0.6))
        trap_2 = ImageMobject("trap_demo_2.png").scale(0.5)
        trap_2.add(Tex("Fig 4b: L=3").next_to(trap_2, DOWN, buff=SMALL_BUFF).scale(0.6))
        trap_3 = ImageMobject("trap_demo_3.png").scale(0.5)
        trap_3.add(Tex("Fig 4c: L=4").next_to(trap_3, DOWN, buff=SMALL_BUFF).scale(0.6))
        trap_plots = Group(trap_1, trap_2, trap_3).arrange(RIGHT)
        trap_text = Tex("Fig 4: Escape time plots with increasing L. (E=1.1)").next_to(trap_plots, DOWN).scale(0.8)
        escape_time_group = Group(trap_plots, trap_text).to_edge(UP)
        self.play(FadeIn(escape_time_group))
        self.wait(0.2)
        self.next_section()
        text = Tex("A lesser amount of trajectories\\ seem to escape for higher L").to_edge(DOWN, buff=2*LARGE_BUFF)
        self.play(Write(text, run_time=0.5))
        self.wait(0.2)
        self.next_section()
        self.play(Unwrite(text, run_time=0.5))
        trap_pass_1 = ImageMobject("trap_pass_demo_1.png").scale(0.5)
        trap_pass_1.add(Tex("Fig 5a: L=2").next_to(trap_pass_1, DOWN, buff=SMALL_BUFF).scale(0.6))
        trap_pass_2 = ImageMobject("trap_pass_demo_2.png").scale(0.5)
        trap_pass_2.add(Tex("Fig 5b: L=3").next_to(trap_pass_2, DOWN, buff=SMALL_BUFF).scale(0.6))
        trap_pass_3 = ImageMobject("trap_pass_demo_3.png").scale(0.5)
        trap_pass_3.add(Tex("Fig 5c: L=4").next_to(trap_pass_3, DOWN, buff=SMALL_BUFF).scale(0.6))
        trap_pass_plots = Group(trap_pass_1, trap_pass_2, trap_pass_3).arrange(RIGHT)
        trap_pass_text = Tex("Fig 5: Escape pass plots with increasing L. (E=1.1)").next_to(trap_pass_plots, DOWN).scale(0.8)
        escape_pass_group = Group(trap_pass_plots, trap_pass_text).to_edge(DOWN)
        self.play(FadeIn(escape_pass_group))
        self.wait(0.2)


class ChaosDefinitionScene(Scene):
    def construct(self):
        slide_title = Tex("Chaotic Dynamics").scale(3)
        what = VGroup(
            Tex("How do we qualify a"),
            Tex("trajectory as chaotic?")
        ).arrange(DOWN).scale(2)
        title_group = VGroup(slide_title, what).arrange(DOWN)
        self.play(Write(title_group[0], run_time=0.5))
        self.wait(0.2)
        self.next_section()
        self.play(Write(title_group[1], run_time=0.5))
        self.wait(0.2)
        self.next_section()
        self.play(Unwrite(title_group, run_time=0.5))
        sensitivity = Tex("Sensitivity to initial conditions").to_edge(UP).scale(1.5)
        self.play(Write(sensitivity, run_time=0.5))
        self.wait(0.2)
        self.next_section()
        stretch = 2
        L = ValueTracker(1.0)
        E = ValueTracker(0.5)
        ax_chaotic = Axes(
            x_length=3 * stretch,
            x_range=(0, 3),
            y_range=(-3, 3)
        ).scale(0.9)
        ax_label_chaotic = ax_chaotic.get_axis_labels(x_label=MathTex('r'), y_label=MathTex('z'))
        r_chaotic_value = 1.25
        r_chaotic = ValueTracker(r_chaotic_value)
        T = ValueTracker(50)
        trajectory_chaotic = MNCTrajectory2D(L, E, r_chaotic, ax_chaotic, T)
        zvc_chaotic = MNCZeroVelocityCurve(L, E, ax_chaotic)
        label_chaotic = Tex("Chaotic Trajectory").next_to(ax_chaotic, DOWN)
        trajectory_plot_chaotic = VGroup(
            ax_chaotic, ax_label_chaotic, label_chaotic, trajectory_chaotic, zvc_chaotic
        ).to_corner(DL)
        self.play(FadeIn(trajectory_plot_chaotic[:-2]), Create(trajectory_plot_chaotic[-1]))
        self.play(Create(trajectory_plot_chaotic[-2]))
        self.wait(0.2)
        self.next_section()
        self.play(r_chaotic.animate.set_value(1.2))
        self.play(r_chaotic.animate.set_value(1.3))
        self.play(r_chaotic.animate.set_value(1.2))
        self.play(r_chaotic.animate.set_value(1.25))
        self.wait(0.2)
        self.next_section()
        ax_ordered = Axes(
            x_length=3 * stretch,
            x_range=(0, 3),
            y_range=(-3, 3)
        ).scale(0.9)
        ax_label_ordered = ax_ordered.get_axis_labels(x_label=MathTex('r'), y_label=MathTex('z'))
        r_ordered_value = 1.05
        r_ordered = ValueTracker(r_ordered_value)
        T = ValueTracker(50)
        trajectory_ordered = MNCTrajectory2D(L, E, r_ordered, ax_ordered, T)
        zvc_ordered = MNCZeroVelocityCurve(L, E, ax_ordered)
        label_ordered = Tex("Ordered Trajectory").next_to(ax_ordered, DOWN)
        trajectory_plot_ordered = VGroup(
            ax_ordered, ax_label_ordered, label_ordered, trajectory_ordered, zvc_ordered
        ).to_corner(DR)
        self.play(FadeIn(trajectory_plot_ordered[:-2]), Create(trajectory_plot_ordered[-1]))
        self.play(Create(trajectory_plot_ordered[-2]))
        self.wait(0.2)
        self.next_section()
        self.play(r_ordered.animate.set_value(0.9))
        self.play(r_ordered.animate.set_value(1.05))
        self.play(r_ordered.animate.set_value(0.9))
        self.play(r_ordered.animate.set_value(1.0))
        self.wait(0.2)


class SaliScene(Scene):
    def construct(self):
        how = Tex("How do we quantify chaos?").scale(2)
        answer = Tex("Single Alignment Index (SALI)").scale(2)
        group = VGroup(how, answer).arrange(DOWN, buff=LARGE_BUFF)
        self.play(Write(group[0], run_time=0.5))
        self.wait(0.2)
        self.next_section()
        self.play(Write(group[1], run_time=0.5))
        self.wait(0.2)
        self.next_section()
        self.play(FadeOut(how, shift=UP), answer.animate.to_edge(UP, buff=SMALL_BUFF).scale(2/3))
        self.wait(0.2)
        self.next_section()
        a = 0.3
        v = 2
        ordered_particle_1 = Dot(color=BLUE).move_to(5.4 * LEFT + 0.6 * UP)
        ordered_particle_2 = Dot(color=BLUE).move_to(6 * LEFT)
        ordered_particle_3 = Dot(color=BLUE).move_to(5.4 * LEFT + 0.6 * DOWN)
        chaotic_particle_1 = Dot(color=BLUE).move_to((6 - 6 * a * np.exp(-a * 6)) * LEFT + 6 * a * np.exp(-a * 6) * UP)
        chaotic_particle_2 = Dot(color=BLUE).move_to(6 * LEFT)
        chaotic_particle_3 = Dot(color=BLUE).move_to((6 - 6 * a * np.exp(-a * 6)) * LEFT + 6 * a * np.exp(-a * 6) * DOWN)

        ordered_deviation_1 = Line(
            stroke_color=WHITE, stroke_width=10
        ).put_start_and_end_on(
            ordered_particle_1.get_center(), ordered_particle_2.get_center()
        ).set_z_index(ordered_particle_1.z_index - 1)
        ordered_deviation_2 = Line(
            stroke_color=WHITE, stroke_width=10
        ).put_start_and_end_on(
            ordered_particle_3.get_center(), ordered_particle_2.get_center()
        ).set_z_index(ordered_particle_1.z_index - 1)

        ordered_deviation_1.add_updater(
            lambda mob: mob.put_start_and_end_on(ordered_particle_1.get_center(), ordered_particle_2.get_center())
        )
        ordered_deviation_2.add_updater(
            lambda mob: mob.put_start_and_end_on(ordered_particle_3.get_center(), ordered_particle_2.get_center())
        )
        ordered_text = Tex("Ordered").next_to(ordered_particle_2, LEFT)
        ordered_group = VGroup(
            ordered_text,
            ordered_particle_1,
            ordered_particle_2,
            ordered_particle_3,
            TracedPath(ordered_particle_1.get_center, stroke_color=ordered_particle_1.get_color(), stroke_width=10),
            TracedPath(ordered_particle_2.get_center, stroke_color=ordered_particle_2.get_color(), stroke_width=10),
            TracedPath(ordered_particle_3.get_center, stroke_color=ordered_particle_3.get_color(), stroke_width=10),
            ordered_deviation_1,
            ordered_deviation_2
        ).to_corner(UL).shift(DOWN)

        chaotic_deviation_1 = Line(
            stroke_color=WHITE, stroke_width=10
        ).put_start_and_end_on(
            chaotic_particle_1.get_center(), chaotic_particle_2.get_center()
        ).set_z_index(chaotic_particle_1.z_index - 1)
        chaotic_deviation_2 = Line(
            stroke_color=WHITE, stroke_width=10
        ).put_start_and_end_on(
            chaotic_particle_3.get_center(), chaotic_particle_2.get_center()
        ).set_z_index(chaotic_particle_1.z_index - 1)

        chaotic_deviation_1.add_updater(
            lambda mob: mob.put_start_and_end_on(chaotic_particle_1.get_center(), chaotic_particle_2.get_center())
        )

        chaotic_deviation_2.add_updater(
            lambda mob: mob.put_start_and_end_on(chaotic_particle_3.get_center(), chaotic_particle_2.get_center())
        )

        chaotic_text = Tex("Chaotic").next_to(chaotic_particle_2, LEFT)
        chaotic_group = VGroup(
            chaotic_text,
            chaotic_particle_1,
            chaotic_particle_2,
            chaotic_particle_3,
            TracedPath(chaotic_particle_1.get_center, stroke_color=chaotic_particle_1.get_color(), stroke_width=10),
            TracedPath(chaotic_particle_2.get_center, stroke_color=chaotic_particle_2.get_color(), stroke_width=10),
            TracedPath(chaotic_particle_3.get_center, stroke_color=chaotic_particle_3.get_color(), stroke_width=10),
            chaotic_deviation_1,
            chaotic_deviation_2
        ).to_corner(DL).shift(UP)
        self.add(ordered_group, chaotic_group)
        self.wait(0.2)
        self.next_section()

        def ordered_updater(mob, dt):
            dx = v * dt
            mob.shift(RIGHT * dx)

        ordered_particle_1.add_updater(ordered_updater)
        ordered_particle_2.add_updater(ordered_updater)
        ordered_particle_3.add_updater(ordered_updater)

        def chaotic_updater_1(mob, dt):
            x, y, z = mob.get_center()
            dx = v * dt
            dy = 2 * v * a ** 2 * np.exp(a * x) * dt
            mob.shift(RIGHT * dx + UP * dy)

        def chaotic_updater_3(mob, dt):
            x, y, z = mob.get_center()
            dx = v * dt
            dy = -2 * v * a ** 2 * np.exp(a * x) * dt
            mob.shift(RIGHT * dx + UP * dy)

        chaotic_particle_1.add_updater(chaotic_updater_1)
        chaotic_particle_2.add_updater(ordered_updater)
        chaotic_particle_3.add_updater(chaotic_updater_3)

        self.wait(4)
        self.next_section()
        self.remove(ordered_group, chaotic_group)
        self.add(ordered_text, chaotic_text)
        ordered_text_description = Tex("- deviation vectors tend to stay unaligned over time").next_to(ordered_text, RIGHT)

        chaotic_text_description = Tex("- deviation vectors tend to align over time").next_to(chaotic_text, RIGHT)
        self.play(
            *[Write(i, run_time=0.5) for i in (ordered_text_description, chaotic_text_description)]
        )
        self.wait(0.2)
        self.next_section()
        deviation_vector_text = VGroup(
            Tex("The deviation vector $\\vec{v}$ has 4 components"),
            MathTex("\\vec{v}=(v_r, v_z, v_{pr}, v_{pz})"),
            Tex("$\\vec{v}$ is time-evolved by the linearized EoM"),
            MathTex("\\dot{v}_i=\\sum_j\\frac{\\partial\\dot x_i}{\\partial x_j}\\vec{v}_j"),
            MathTex(r'\dot{p}_r={L\over r^3}-\frac{r}{4}-\frac{r}{R^3}'),
            MathTex("\\dot{v}_{pr}=-\\left(\\frac{1}{4}+\\frac{3L^2}{r^4}+\\frac{1}{R^3}-\\frac{3 r^2}{R^5}\\right)v_r+\\left(\\frac{3rz}{R^5}\\right)v_z")
        ).arrange(DOWN).next_to(answer, DOWN)
        label20 = MathTex("(20)").to_edge(RIGHT)
        label4 = MathTex("(4)").to_edge(RIGHT)
        label21 = MathTex("(21)").to_edge(RIGHT)
        label22 = MathTex("(22)").to_edge(RIGHT)
        label20.move_to(label20.get_x()*RIGHT+deviation_vector_text[1].get_y()*UP)
        deviation_vector_text[1].add(label20)
        label21.move_to(label20.get_x()*RIGHT+deviation_vector_text[3].get_y()*UP)
        deviation_vector_text[3].add(label21)
        label4.move_to(label20.get_x()*RIGHT+deviation_vector_text[4].get_y()*UP)
        deviation_vector_text[4].add(label4)
        label22.move_to(label20.get_x()*RIGHT+deviation_vector_text[5].get_y()*UP)
        deviation_vector_text[5].add(label22)
        self.play(
            Unwrite(VGroup(ordered_text, ordered_text_description, chaotic_text, chaotic_text_description), run_time=0.5),
            Write(deviation_vector_text[:2], run_time=0.5)
        )
        self.wait(0.2)
        self.next_section()
        for i in deviation_vector_text[2:]:
            self.play(
                Write(i, run_time=0.5)
            )
            self.wait(0.2)
            self.next_section()
        self.play(Unwrite(deviation_vector_text, run_time=0.5))
        sali_text = VGroup(
            Tex("SALI can be calculated using 2 arbitrary deviation vectors"),
            MathTex("SALI(t)=\\log_{10}(\\min(||\\hat{v}_1+\\hat{v}_2||, ||\\hat{v}_1-\\hat{v}_2||))"),
            Tex("For chaotic trajectories,"),
            MathTex("SALI(t\\to\infty)\\to-\\infty"),
            Tex("For ordered trajectories,"),
            MathTex("SALI(t\\to\infty)\\to\\text{finite}")
        ).arrange(DOWN)
        label23 = MathTex("(23)").to_edge(RIGHT)
        label24 = MathTex("(24)").to_edge(RIGHT)
        label25 = MathTex("(25)").to_edge(RIGHT)
        label23.move_to(label23.get_x()*RIGHT+sali_text[1].get_y()*UP)
        sali_text[1].add(label23)
        label24.move_to(label23.get_x()*RIGHT+sali_text[3].get_y()*UP)
        sali_text[3].add(label24)
        label25.move_to(label23.get_x()*RIGHT+sali_text[5].get_y()*UP)
        sali_text[5].add(label25)
        self.play(Write(sali_text[:2], run_time=0.5))
        self.wait(0.2)
        self.next_section()
        for i in sali_text[2:]:
            self.play(
                Write(i, run_time=0.5)
            )
            self.wait(0.2)
            self.next_section()
        self.play(Unwrite(sali_text, run_time=0.5), Unwrite(answer, run_time=0.5))
        self.wait(0.2)


class SaliPlot(Scene):
    def construct(self):
        SALI_plot = ImageMobject("sali_plot.png").scale(0.8)
        SALI_plot_scale = ImageMobject("sali_plot_scale.png").scale(0.8)
        SALI_plot_text = Tex(r"Fig 6: SALI plot\\ (L, E)=(1, 0.42)").next_to(SALI_plot, DOWN)
        SALI_axis = Axes(
            x_range=(0.5177782924154619, 1.7573551332081243),
            y_range=(-0.8609681374459811, 0.860968137445981),
            x_length=SALI_plot_scale.width,
            y_length=SALI_plot_scale.height
        ).shift(
            0.06 * UP * SALI_plot_scale.height +
            0.0185 * LEFT * SALI_plot_scale.width
        )
        SALI_group = Group(SALI_axis, SALI_plot, SALI_plot_text).to_edge(LEFT)
        self.play(FadeIn(SALI_group[1:]))
        self.wait(0.2)
        self.next_section()
        stretch = 2
        traj_axis = Axes(
            x_length=3 * stretch,
            x_range=(0, 2),
            y_range=(-3, 3)
        )
        ax_label = traj_axis.get_axis_labels(x_label=MathTex('r'), y_label=MathTex('z'))
        L = ValueTracker(1.0)
        E = ValueTracker(0.42)
        r = ValueTracker(1.0)
        pr = ValueTracker(0)
        T = ValueTracker(50)
        zvc = MNCZeroVelocityCurve(L, E, traj_axis)
        trajectory = MNCTrajectory2D(L, E, r, traj_axis, T, pr)
        trajectory_plot = VGroup(traj_axis, ax_label, trajectory, zvc).to_edge(RIGHT)
        ic_point = VGroup(
            Dot(SALI_axis.c2p(r.get_value(), pr.get_value()), stroke_width=3, color=BLACK, stroke_color=WHITE),
            Circle(radius=0.5, stroke_opacity=0.5, color=BLACK).move_to(SALI_axis.c2p(r.get_value(), pr.get_value()))
        )
        ic_point.add_updater(
            lambda mob: mob.move_to(SALI_axis.c2p(r.get_value(), pr.get_value()))
        )
        self.play(Create(trajectory_plot), Create(ic_point))
        self.wait(0.2)
        self.next_section()
        self.play(r.animate.set_value(0.98))
        self.wait(0.5)
        self.play(r.animate.set_value(1.15))
        self.wait(0.5)
        self.play(r.animate.set_value(0.75), pr.animate.set_value(-0.4))
        self.wait(0.2)
        self.next_section()
        self.play(r.animate.set_value(0.7), pr.animate.set_value(-0.65))
        self.wait(0.5)
        self.play(r.animate.set_value(1.1), pr.animate.set_value(-0.72))
        self.wait(0.5)
        self.play(r.animate.set_value(1), pr.animate.set_value(-0.34))
        self.wait(0.5)
        self.play(r.animate.set_value(0.95), pr.animate.set_value(-0.24))
        self.wait(0.5)
        self.play(r.animate.set_value(0.77), pr.animate.set_value(0))
        self.wait(0.2)


class SaliTrap(Scene):
    def construct(self):
        trap_1 = ImageMobject("trap_demo_1.png").scale(0.5)
        trap_1.add(Tex("Fig 4a: L=2").next_to(trap_1, DOWN, buff=SMALL_BUFF).scale(0.6))
        trap_2 = ImageMobject("trap_demo_2.png").scale(0.5)
        trap_2.add(Tex("Fig 4b: L=3").next_to(trap_2, DOWN, buff=SMALL_BUFF).scale(0.6))
        trap_3 = ImageMobject("trap_demo_3.png").scale(0.5)
        trap_3.add(Tex("Fig 4c: L=4").next_to(trap_3, DOWN, buff=SMALL_BUFF).scale(0.6))
        trap_plots = Group(trap_1, trap_2, trap_3).arrange(RIGHT)
        trap_text = Tex("Fig 4: Escape time plots with increasing L. (E=1.1)").next_to(trap_plots, DOWN).scale(0.8)
        trap_plots = Group(trap_plots, trap_text).to_edge(DOWN)

        sali_1 = ImageMobject("sali_trap_1.png").scale(0.5)
        sali_1.add(Tex("Fig 7a: L=2").next_to(sali_1, DOWN, buff=SMALL_BUFF).scale(0.6))
        sali_2 = ImageMobject("sali_trap_2.png").scale(0.5)
        sali_2.add(Tex("Fig 7b: L=3").next_to(sali_2, DOWN, buff=SMALL_BUFF).scale(0.6))
        sali_3 = ImageMobject("sali_trap_3.png").scale(0.5)
        sali_3.add(Tex("Fig 7c: L=4").next_to(sali_3, DOWN, buff=SMALL_BUFF).scale(0.6))
        sali_plots = Group(sali_1, sali_2, sali_3).arrange(RIGHT)
        sali_text = Tex("Fig 7: SALI plots with increasing L. (E=1.1)").next_to(sali_plots, DOWN).scale(0.8)
        sali_group = Group(sali_plots, sali_text).to_edge(UP)

        self.play(FadeIn(sali_group))
        self.wait(0.2)
        self.next_section()
        self.play(FadeIn(trap_plots))
        self.wait(0.2)


class SaliTrapPlot(Scene):
    def construct(self):
        SALI_plot = ImageMobject("sali_trap_3.png")
        SALI_plot_scale = ImageMobject("sali_trap_3_scale.png")
        SALI_plot_text = Tex(r"Fig 7c: SALI plot\\ (L, E)=(4, 1.1)").next_to(SALI_plot, DOWN)
        SALI_axis = Axes(
            x_range=(1.9892831554621069, 3.698518202727969),
            y_range=(-0.8744109279194485, 0.8744109279194486),
            x_length=SALI_plot_scale.width,
            y_length=SALI_plot_scale.height
        ).shift(
            0.09 * UP * SALI_plot_scale.height +
            0.008 * LEFT * SALI_plot_scale.width
        ).set_color(RED)
        SALI_group = Group(SALI_axis, SALI_plot, SALI_plot_text).to_edge(LEFT)
        self.play(FadeIn(SALI_group[1:]))
        stretch = 2
        traj_axis = Axes(
            x_length=3 * stretch,
            x_range=(1.8, 3.8),
            y_range=(-10, 10, 5)
        )
        ax_label = traj_axis.get_axis_labels(x_label=MathTex('r'), y_label=MathTex('z'))
        L = ValueTracker(4.0)
        E = ValueTracker(1.1)
        r = ValueTracker(3.3)
        pr = ValueTracker(0.)
        T = ValueTracker(300)
        zvc = MNCZeroVelocityCurve(L, E, traj_axis)
        trajectory = MNCTrajectory2D(L, E, r, traj_axis, T, pr)
        asymptotes = MNCAsymptotes(L, E, traj_axis)
        trajectory_plot = VGroup(traj_axis, ax_label, trajectory, asymptotes, zvc).to_edge(RIGHT)
        ic_point = VGroup(
            Dot(SALI_axis.c2p(r.get_value(), pr.get_value()), stroke_width=3, color=BLACK, stroke_color=WHITE),
            Circle(radius=0.5, stroke_opacity=0.5, color=BLACK).move_to(SALI_axis.c2p(r.get_value(), pr.get_value()))
        )
        ic_point.add_updater(
            lambda mob: mob.move_to(SALI_axis.c2p(r.get_value(), pr.get_value()))
        )
        self.play(Create(trajectory_plot), Create(ic_point))
        self.wait(0.2)
        self.next_section()
        self.play(r.animate.set_value(3.5))
        self.wait(0.2)
        self.next_section()
        self.play(r.animate.set_value(3.6))
        self.wait(0.2)
        self.next_section()
        self.play(r.animate.set_value(3.0))
        self.wait(0.2)
        self.next_section()
        self.play(r.animate.set_value(2.8))
        self.wait(0.2)


class CountOrdered(Scene):
    def construct(self):
        text = VGroup(
            Tex("For what (L, E) values does these \"trapped\" orbits exist in?"),
            Tex("Set SALI$>$-1 as ordered, and SALI$<$-1 as chaotic"),
            Tex("Proportion of ordered trajectories\\ for a given (L, E)")
        ).arrange(DOWN).to_edge(UP)

        order_plot = ImageMobject("order_plot.png").scale(0.9)
        order_plot_text = Tex(r"Fig 8: \% of ordered trajectories given (L, E)").next_to(order_plot, DOWN)
        order_group = Group(order_plot, order_plot_text).to_edge(DOWN)
        for i in text:
            self.play(Write(i, run_time=0.5))
            self.wait(0.2)
            self.next_section()
        self.play(FadeIn(order_group))
        self.wait(0.2)
        self.next_section()
        self.play(FadeOut(order_group, run_time=0.5), FadeOut(text, run_time=0.5))
        self.wait(0.2)
        self.next_section()
        mnc_qualities = Tex("MNC qualities").scale(2).to_edge(UP)
        chaotic_list = BulletedList(
            r'There exists "trapped orbits",\\ordered trajectories with\\escaping energies that\\will never escape',
            r"More orbits are trapped\\as the magnitude of\\L increases"
        ).scale(0.9)
        escape_list = BulletedList(
            r"For E greater than 1, trajectories\\can escape along $z\pm\infty$",
            r"There exists initial conditions\\with escapes that are\\insensitive to perturbations",
            r"There exists initial conditions\\with escapes that are\\sensitive to perturbations",
            r"Initial conditions on boundaries\\of the two aformentioned regions\\have high escape times"
        ).scale(0.9).to_edge(RIGHT)
        self.play(Write(chaotic_list, run_time=0.5), Write(mnc_qualities, run_time=0.5))
        self.wait(0.2)
        self.next_section()
        self.play(
            chaotic_list.animate.to_edge(LEFT),
            mnc_qualities.animate.to_corner(UL)
        )
        self.play(Write(escape_list, run_time=0.5))
        self.wait(0.2)


class MBHamiltonianScene(Scene):
    def construct(self):
        TitleScreen = Tex(r"Magnetized Black\\ Hole").scale(2.5)
        self.play(Write(TitleScreen, run_time=0.5))
        self.wait(0.2)
        self.next_section()
        super_ham = MathTex(
            r"\mathcal{H}=\frac{1}{2m}\left(g^{\mu\nu}(p_\mu-qA_\mu)(p_\nu-qA_\nu)+(mc)^2\right)"
        )
        metric = MathTex(
            r"g_{\mu\nu}=\text{diag}\left(-fc^2, {1\over f}, R^2, R^2\sin^2\theta\right)"
        )
        f = MathTex("f=1-{2GM \over Rc^2}")
        A = MathTex(
            r"A^\mu=\left(0,0,0, {B\over 2}\right)"
        )
        group_1 = VGroup(super_ham, metric, f, A).arrange(DOWN)
        label26 = MathTex('(26)').next_to(super_ham, RIGHT).to_edge(RIGHT)
        label27 = MathTex('(27)').next_to(metric, RIGHT).to_edge(RIGHT)
        label28 = MathTex('(28)').next_to(f, RIGHT).to_edge(RIGHT)
        label29 = MathTex('(29)').next_to(A, RIGHT).to_edge(RIGHT)
        self.play(Unwrite(TitleScreen, run_time=0.5))
        self.play(*[Write(i, run_time=0.5) for i in (super_ham, label26)])
        self.wait(0.2)
        self.next_section()
        self.play(*[Write(i, run_time=0.5) for i in (metric, label27)])
        self.play(*[Write(i, run_time=0.5) for i in (f, label28)])
        self.wait(0.2)
        self.next_section()
        self.play(*[Write(i, run_time=0.5) for i in (A, label29)])
        self.wait(0.2)
        self.next_section()
        self.play(*[Unwrite(i, run_time=0.5) for i in (group_1, label26, label27, label28, label29)])
        ham = MathTex(
            r"\mathcal{H}=\frac{1}{2}\left(-{", r"\mathcal{E}", r"^2\over f}+fp_R^2+{p_\theta^2\over R^2}+\left(", r"{L\over R\sin\theta}-",
            "B", r"R\sin\theta\right)^2+1\right)"
        ).scale(0.9).to_corner(UL).shift(2*DOWN)
        label30 = MathTex('(30)').next_to(ham, RIGHT).to_edge(RIGHT)
        char_l = MathTex(r"l_c=2GM/c^2").to_corner(DL).shift(RIGHT * 1.8 + 2.5 * UP)
        label31 = MathTex('(31)').next_to(char_l, RIGHT, buff=1)
        char_t = MathTex(r"t_c=l_c/c").next_to(label31, RIGHT, buff=2)
        label32 = MathTex('(32)').next_to(char_t, RIGHT, buff=1)
        self.play(*[Write(i) for i in (ham, label30, f)], run_time=0.5)
        self.wait(0.2)
        self.next_section()
        self.play(*[Write(i) for i in (char_t, char_l, label31, label32)], run_time=1)
        self.wait(0.2)
        self.next_section()
        ham[1].set_color(GREEN)
        self.play(Indicate(ham[1], color=GREEN, scale_factor=1.5), run_time=1.5)
        self.wait(0.2)
        self.next_section()
        ham[3][0].set_color(GREEN)
        self.play(Indicate(ham[3][0], color=GREEN, scale_factor=1.5), run_time=1.5)
        self.wait(0.2)
        self.next_section()
        ham[4].set_color(GREEN)
        self.play(Indicate(ham[4], color=GREEN, scale_factor=1.5), run_time=1.5)
        self.wait(0.2)
        self.next_section()
        all_text = VGroup()
        for i in self.mobjects:
            if isinstance(i, (MathTex, VGroup)):
                all_text.add(i)
        self.play(Unwrite(all_text, run_time=0.5))

        normalization = MathTex(r"p^\mu p_\mu=-1")
        energy_eqn = MathTex(
            r"\mathcal{E}^2=f^2p_R^2+{fp_\theta\over R^2}+U"
        )
        potential = MathTex(r"U=\left(1-{1\over R}\right)\left(1+\left({L\over r}-Br\right)^2\right)")
        group = VGroup(normalization, energy_eqn, potential).arrange(DOWN)
        normalization.add(MathTex('(33)').next_to(normalization, RIGHT).to_edge(RIGHT))
        energy_eqn.add(MathTex('(34)').next_to(energy_eqn, RIGHT).to_edge(RIGHT))
        potential.add(MathTex('(35)').next_to(potential, RIGHT).to_edge(RIGHT))
        for i in group:
            self.play(Write(i, run_time=0.5))
            self.wait(0.2)
            self.next_section()
        self.play(Unwrite(group, run_time=0.5))
        self.wait(0.2)


class MBHZVC(Scene):
    def construct(self):
        stretch = 2
        ax = Axes(
            x_length=3 * stretch,
            x_range=(0, 6.5)
        )
        ax_label = ax.get_axis_labels(x_label=MathTex('r'), y_label=MathTex('z'))
        r_esc = ValueTracker(3.7)
        p = ValueTracker(1.8)
        E = ValueTracker(0.9)
        zvc = MBHZeroVelocityCurve(r_esc, p, E, ax)
        asymptotes = MBHAsymptotes(r_esc, p, E, ax)

        text_B = DecimalNumber(
            MBHreparametrization(r_esc.get_value(), p.get_value(), E.get_value())[0],
            num_decimal_places=3
        ).move_to(4 * LEFT + 1 * UP)
        text_L = DecimalNumber(
            MBHreparametrization(r_esc.get_value(), p.get_value(), E.get_value())[1],
            num_decimal_places=3
        ).next_to(text_B, DOWN)
        text_E = DecimalNumber(
            MBHreparametrization(r_esc.get_value(), p.get_value(), E.get_value())[2],
            num_decimal_places=3
        ).next_to(text_L, DOWN)
        text_B.add_updater(lambda mob: mob.set_value(MBHreparametrization(r_esc.get_value(), p.get_value(), E.get_value())[0]))
        text_L.add_updater(lambda mob: mob.set_value(MBHreparametrization(r_esc.get_value(), p.get_value(), E.get_value())[1]))
        text_E.add_updater(lambda mob: mob.set_value(MBHreparametrization(r_esc.get_value(), p.get_value(), E.get_value())[2]))
        B_label = MathTex("B=").next_to(text_B, LEFT)
        L_label = MathTex("L=").next_to(text_L, LEFT)
        E_label = MathTex(r"\mathcal{E}=").next_to(text_E, LEFT)
        BHcurve = ParametricFunction(
            lambda t: (ax.c2p(np.cos(t), 0)[0], ax.c2p(0, np.sin(t))[1], 0), t_range=(-PI/2, PI/2), color=ORANGE
        )
        self.add(ax, ax_label, zvc, asymptotes, text_B, text_L, text_E, B_label, L_label, E_label, BHcurve)

        self.wait(0.2)
        self.next_section()
        self.play(E.animate.set_value(zvc.values.get_U_eq()[1]*1.001))
        self.wait(0.2)
        self.next_section()
        self.play(E.animate.set_value(1.01))
        self.wait(0.2)
        self.next_section()
        self.play(E.animate.set_value(1.3))
        self.wait(0.2)
        self.next_section()
        self.play(E.animate.set_value(1.01))
        self.wait(0.2)
        self.next_section()
        self.play(E.animate.set_value(0.75))
        self.play(r_esc.animate.set_value(2.0))
        self.wait(0.2)


class MBHRescale(Scene):
    def construct(self):
        reparam = Tex("Reparametrization").scale(1.5).to_edge(UP)
        E0 = MathTex(r"(B, L, \mathcal{E})\Rightarrow (r_{esc}, \rho, E)")
        E1 = MathTex(r"\mathcal{E}^2=U_{esc,min}E=\left[1+2B\left(|L|-L\right)\right]E")
        E2 = MathTex(r"L_\pm=\pm\frac{r_{\text{esc}}\sqrt[4]{(3-\rho)(3\rho-1)}}{\sqrt{2\left(4\rho^2-9\rho+3\pm\sqrt{(3\rho-1)(3-\rho)}\right)}}")
        E3 = MathTex(r"B_\pm=\frac{\sqrt[4]{(3-\rho)(3\rho-1)}}{r_{\text{esc}}\sqrt{2\left(4\rho^2-9\rho+3\pm\sqrt{(3\rho-1)(3-\rho)}\right)}}")
        text_group = VGroup(reparam, E0, E1, E2, E3).arrange(DOWN)
        label28 = MathTex('(36)').next_to(E0, RIGHT).to_edge(RIGHT)
        label29 = MathTex('(37)').next_to(E1, RIGHT).to_edge(RIGHT)
        label30 = MathTex('(38)').next_to(E2, RIGHT).to_edge(RIGHT)
        label31 = MathTex('(39)').next_to(E3, RIGHT).to_edge(RIGHT)
        E0.add(label28)
        E1.add(label29)
        E2.add(label30)
        E3.add(label31)
        self.play(Write(text_group[0], run_time=0.5))
        self.wait(0.2)
        self.next_section()
        self.play(Write(text_group[1], run_time=0.5))
        self.wait(0.2)
        self.next_section()
        self.play(Write(text_group[2:], run_time=0.5))
        self.wait(0.2)
        self.next_section()
        self.play(Unwrite(text_group[1:], run_time=0.5))
        E4 = Tex(r"For a given ($r_{esc}, \rho, E, R_0, p_{R0}, \theta_0$),")
        E5 = Tex("trajectories are symmetrical under")
        E6 = MathTex(r"(B_+, L_+, \mathcal{E})\Leftrightarrow (B_-, L_-, \mathcal{E}).")
        E7 = Tex(r"All possible values of ($B, L, \mathcal{E}$)")
        E8 = Tex(r"are included in ($r_{esc}, \rho, E$)")
        text_group_2 = VGroup(E4, E5, E6, E7, E8).arrange(DOWN, buff=MED_LARGE_BUFF)
        label32 = MathTex('(40)').next_to(E6, RIGHT).to_edge(RIGHT)
        E6.add(label32)
        self.play(Write(text_group_2[:2], run_time=0.5))
        self.wait(0.2)
        self.next_section()
        self.play(Write(text_group_2[2], run_time=0.5))
        self.wait(0.2)
        self.next_section()
        self.play(Write(text_group_2[3:], run_time=0.5))
        self.wait(0.2)


class MBHEscapePlot(Scene):
    def construct(self):
        title = Tex("Escape\\ Dynamics").scale(3)
        self.play(Write(title, run_time=0.5))
        self.wait(0.2)
        self.next_section()
        self.play(Unwrite(title, run_time=0.5))
        self.wait(0.2)
        self.next_section()
        escape_time = ImageMobject("MBH_time.png").scale(0.8)
        escape_time_scale = ImageMobject("MBH_time_scale.png").scale(0.8)
        escape_time_text = Tex(r"Fig 9: Escape time plot\\ $(r_{esc}, \rho, E)=(1.8, 1.2, 1.2)$").next_to(escape_time, DOWN)
        escape_axis = Axes(
            x_range=(1.0038801161495363, 2.164034844860891),
            y_range=(-1.088166409927207, 1.088166409927207),
            x_length=escape_time_scale.width,
            y_length=escape_time_scale.height
        ).shift(
            0.09 * UP * escape_time_scale.height +
            0.015 * LEFT * escape_time_scale.width
        ).set_color(RED)
        escape_time_group = Group(escape_time, escape_time_text, escape_axis).to_edge(LEFT)
        self.play(FadeIn(escape_time_group[:2]))
        self.wait(0.2)
        self.next_section()
        stretch = 2
        traj_axis = Axes(
            x_length=3 * stretch,
            x_range=(0, 2.5),
            y_range=(-5, 5)
        )
        ax_label = traj_axis.get_axis_labels(x_label=MathTex('r'), y_label=MathTex('z'))
        r_esc = ValueTracker(1.8)
        p = ValueTracker(1.2)
        E = ValueTracker(1.2)
        T = ValueTracker(100)
        r = ValueTracker(2)
        pr = ValueTracker(0.)
        zvc = MBHZeroVelocityCurve(r_esc, p, E, traj_axis)
        asymptotes = MBHAsymptotes(r_esc, p, E, traj_axis)
        trajectory = MBHTrajectory2D(r_esc, p, E, r, traj_axis, T, pr)
        trajectory_plot = VGroup(traj_axis, ax_label, trajectory, zvc, asymptotes).to_edge(RIGHT)
        ic_point = VGroup(
            Dot(escape_axis.c2p(r.get_value(), pr.get_value()), color=BLACK),
            Circle(radius=0.5, stroke_opacity=0.5, color=BLACK).move_to(escape_axis.c2p(r.get_value(), pr.get_value()))
        )
        ic_point.add_updater(
            lambda mob: mob.move_to(escape_axis.c2p(r.get_value(), pr.get_value()))
        )
        self.play(Create(trajectory_plot), Create(ic_point))
        self.wait(0.2)
        self.next_section()
        self.play(r.animate.set_value(1.66), pr.animate.set_value(0.2))
        self.wait(0.2)
        self.next_section()
        self.play(pr.animate.set_value(0.5))
        self.wait(0.2)
        self.next_section()
        self.play(pr.animate.set_value(0.476), run_time=3)
        self.wait(0.2)
        self.next_section()
        new_axis = Axes(
                x_length=3 * stretch,
                x_range=(0, 2.5),
                y_range=(-0.2, 0.2)
            ).to_edge(RIGHT)
        new_label = new_axis.get_axis_labels(x_label=MathTex('r'), y_label=MathTex('z'))
        new_zvc = MBHZeroVelocityCurve(r_esc, p, E, new_axis)
        new_asymptotes = MBHAsymptotes(r_esc, p, E, new_axis)
        new_trajectory = MBHTrajectory2D(r_esc, p, E, r, new_axis, T, pr)
        self.play(r.animate.set_value(1.2655), pr.animate.set_value(0.056))
        self.wait(0.2)
        new_trajectory_plot = VGroup(new_axis, new_label, new_trajectory, new_zvc, new_asymptotes).to_edge(RIGHT)
        self.remove(trajectory_plot)
        T.set_value(10)
        self.add(new_trajectory_plot)
        self.wait(0.2)
        self.next_section()

        escape_time_zoomed = ImageMobject("MBH_time_zoomed.png").scale(0.8)
        escape_time_zoomed_scale = ImageMobject("MBH_time_zoomed_scale.png").scale(0.8)
        escape_time_zoomed_text = Tex(r"Fig 10: Zoomed escape time\\ plot $(r_{esc}, \rho, E)=(1.8, 1.2, 1.2)$").next_to(escape_time_zoomed, DOWN)
        escape_axis_zoomed = Axes(
            x_range=(1.265, 1.27),
            y_range=(0.055, 0.06),
            x_length=escape_time_zoomed_scale.width,
            y_length=escape_time_zoomed_scale.height
        ).shift(
            0.09 * UP * escape_time_zoomed_scale.height +
            0.0 * LEFT * escape_time_zoomed_scale.width
        )
        escape_time_zoomed_group = Group(escape_axis_zoomed, escape_time_zoomed, escape_time_zoomed_text).to_edge(LEFT)
        ic_point.clear_updaters().set_z_index(escape_time_zoomed_group.z_index + 1)
        BHcurve = ParametricFunction(
            lambda t: (new_axis.c2p(np.cos(t), 0)[0], new_axis.c2p(0, np.sin(t))[1], 0), t_range=(-PI/2, PI/2), color=ORANGE
        )
        self.play(
            FadeOut(escape_time_group),
            FadeIn(escape_time_zoomed_group),
            FadeIn(BHcurve),
            ic_point.animate.move_to(escape_axis_zoomed.c2p(r.get_value(), pr.get_value()))
        )
        ic_point.add_updater(
            lambda mob: mob.move_to(escape_axis_zoomed.c2p(r.get_value(), pr.get_value()))
        )
        self.wait(0.2)
        self.next_section()
        self.play(r.animate.set_value(1.2665), run_time=3)
        self.wait(0.2)


class MBHChaos(Scene):
    def construct(self):
        title = Tex("Chaotic Dynamics").scale(3)
        self.play(Write(title, run_time=0.5))
        self.wait(0.2)
        self.next_section()
        self.play(Unwrite(title, run_time=0.5))
        trap_1 = ImageMobject("MBH_order_plot_1.png").scale(0.6)
        trap_1.add(Tex("Fig 11a: $\\rho=1.2$").next_to(trap_1, DOWN))
        trap_2 = ImageMobject("MBH_order_plot_2.png").scale(0.6)
        trap_2.add(Tex("Fig 11b: $\\rho=1.4$").next_to(trap_2, DOWN))
        trap_3 = ImageMobject("MBH_order_plot_3.png").scale(0.6)
        trap_3.add(Tex("Fig 11c: $\\rho=1.6$").next_to(trap_3, DOWN))
        trap_plots = Group(trap_1, trap_2, trap_3).arrange(RIGHT)
        trap_text = Tex("Fig 11: \% of ordered trajectories given ($r_{esc}, \\rho, E$)").next_to(trap_plots, DOWN)
        trap_group = Group(trap_plots, trap_text).to_edge(UP, buff=LARGE_BUFF)
        self.play(FadeIn(trap_group))
        self.wait(0.2)
        self.next_section()
        conclusion = Tex(r"More orbits are trapped as the magnitude of\\$r_{esc}$ increases").scale(1.3).to_edge(DOWN)
        self.play(Write(conclusion, run_time=0.5))
        self.wait(0.2)


class MBHConclusion(Scene):
    def construct(self):
        mbh_qualities = Tex("MBH qualities").scale(2).to_edge(UR)
        chaotic_list = BulletedList(
            r'There exists "trapped orbits",\\ordered trajectories with\\escaping energies that\\will never escape',
            r"More orbits are trapped\\as the magnitude of\\$r_{esc}$ increases",
            r"Almost no orbits are trapped\\ when BH capture is possible"
        ).scale(0.9).to_edge(RIGHT).shift(0.5*DOWN)
        escape_list = BulletedList(
            r"For high enough energies,\\trajectories can escape along\\$z\pm\infty$ or captured by the BH",
            r"There exists initial conditions\\with escapes that are\\insensitive to perturbations",
            r"There exists initial conditions\\with escapes that are\\sensitive to perturbations",
            r"Initial conditions on boundaries\\of the two aformentioned regions\\have high escape times except when\\the ordered region is BH capture"
        ).scale(0.8).to_edge(LEFT)
        self.play(Write(mbh_qualities, run_time=0.5))
        self.wait(0.2)
        self.next_section()
        for i in escape_list:
            self.play(Write(i, run_time=0.5))
            self.wait(0.2)
            self.next_section()
        self.play(Write(chaotic_list[0], run_time=0.5))
        self.wait(0.2)
        self.next_section()
        self.play(Write(chaotic_list[1], run_time=0.5))
        self.wait(0.2)
        self.next_section()
        self.play(Write(chaotic_list[2], run_time=0.5))
        self.wait(0.2)


class Conclusion(Scene):
    def construct(self):
        conc = Tex("Conclusion").scale(1.2).to_corner(UL)
        self.play(Write(conc, run_time=0.5))
        self.wait(0.2)
        list = BulletedList(
            r"For high enough energies, MNC and MBH trajectories\\ can escape.",
            r"MBH trajectories can also be captured by the black hole",
            r"Near regions with ordered escape structure,\\ escape times are generally higher",
            r"There exists initial conditions with escaping energies that\\do not escape, which produces trapped orbits",
            r"Energy is not the only determining factor for escape."
        ).shift(0.5*DOWN).to_edge(LEFT)
        for i in list:
            self.play(Write(i, run_time=0.5))
            self.wait(0.2)
            self.next_section()

