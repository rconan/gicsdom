use console::{style, Term};
use crossbeam_channel::{bounded, unbounded, Receiver, Sender};
use gicsdom::astrotools;
use std::thread;
use std::time::Duration;

struct Dyad<T> {
    send: Sender<T>,
    recv: Receiver<T>,
}
struct Session<S, C> {
    server: Dyad<S>,
    client: Dyad<C>,
}
impl<S, C> Session<S, C> {
    fn new(server_: (Sender<S>, Receiver<S>), client_: (Sender<C>, Receiver<C>)) -> Self {
        Session {
            server: Dyad {
                send: server_.0,
                recv: server_.1,
            },
            client: Dyad {
                send: client_.0,
                recv: client_.1,
            },
        }
    }
    fn dual_channel(&self) -> (Sender<C>, Receiver<S>) {
        (self.client.send.clone(), self.server.recv.clone())
    }
    fn main_channel(&self) -> (Sender<S>, Receiver<C>) {
        (self.server.send.clone(), self.client.recv.clone())
    }
}
impl<T> Drop for Dyad<T>  {
    fn drop(&mut self) {
        drop(&self.send);
        drop(&self.recv);
    }
}
impl<S, C> Drop for Session<S, C> {
    fn drop(&mut self) {
        drop(&self.server);
        drop(&self.client);
    }
}


fn main() {
    let term = Term::buffered_stdout();
    term.clear_screen().unwrap();


    let n_probe = 3;
    let probe_sessions: [Session<astrotools::Observation, bool>; 3] = [
        Session::new(bounded(0), bounded(0)),
        Session::new(bounded(0), bounded(0)),
        Session::new(bounded(0), bounded(0)),
    ];

    let mut h = vec![];

    for (i, session) in probe_sessions.iter().enumerate() {
        let term_c = term.clone();

        let (sprobe, rprobe) = session.dual_channel();

        h.push(thread::spawn(move || loop {
            let obs = rprobe.recv().unwrap();

            term_c
                .write_line(&format!("PROBE #{}  {}", i, obs.datetime))
                .unwrap();

            sprobe.send(true).unwrap();
        }));
    }

    let mut obs = astrotools::Observation::new(
        "2020-03-20T00:00:00",
        astrotools::GMT_LAT,
        astrotools::GMT_LONG,
    );


    for _ in 0..11 {
        obs.add_seconds(1f64);
        term.move_cursor_to(0, 0).unwrap();
        term.write_line(&format!("{:=>75}", "",)).unwrap();
        term.write_line(&format!("TELESCOPE {}", obs.datetime))
            .unwrap();

        probe_sessions.iter().for_each(|x| x.main_channel().0.send(obs.clone()).unwrap());
        probe_sessions.iter().map(|x| x.main_channel().1.recv().unwrap()).for_each(drop);

        term.write_line(&format!("{:=>75}", "",)).unwrap();
        term.flush().unwrap();
        thread::sleep(Duration::from_secs(1));
    }

    for session in probe_sessions.iter() {
        drop(session);
    }
    for h_ in h {
        println!("join");
        h_.join().unwrap();
    }
}
