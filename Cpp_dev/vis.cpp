#include <gtkmm.h>
#include <iostream>

using namespace std;

void on_button_click();

int main (int argc, char* argv[]) {
	Gtk::Main kit(argc, argv);

	Gtk::Window window;
	Gtk::Button button("Click me");

	window.set_default_size(640, 480);
	window.set_title("Super-noobie program");
	window.set_position(Gtk::WIN_POS_CENTER);
	window.set_border_width(10);

	button.signal_clicked().connect(
		sigc::ptr_fun(&on_button_click)
		);
	button.show();
	window.add(button);

	Gtk::Main::run(window);

	return 0;
}

void on_button_click() {
	static int i = 0;
	cout << ++i << endl;
}