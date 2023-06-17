// group.rs
extern crate rug;
use rug::Integer;
use rug::rand::RandState;

pub struct Group {
    pub n: Integer,
    pub p: Integer,
    pub q: Integer,
    pub n_sq: Integer,
    pub n_root: Integer,
    pub g: Integer,
    pub mu: Integer,
    pub phi_n: Integer,
    pub delta: Integer,
    // ... other fields
}

impl Group {
    // A constructor method
    pub fn new(bit_len: u64) -> Group {
        let p;
        let q;
        if bit_len == 10 {
            p = Integer::from(107);
            q = Integer::from(83);
        }else if bit_len == 100 {
            p = Integer::from_str_radix("6122421090493547576937037317561418841225758554253106999", 10).unwrap();
            q = Integer::from_str_radix("5846418214406154678836553182979162384198610505601062333", 10).unwrap();
        } else if bit_len == 3072 {
            p = Integer::from_str_radix("3793070626895572467283950491156608172646280619543651052186002962401960948344595436008317725360316523061862837176776937156394891603251329682527611995573870675378784146300080405741531439158689694869658339169636599335667867300511737877416759591946775631328177983594319449174751169843162114768684946303385589552432358881413822213266477101325470478323250412128032770624674451411359365977206532598270731108378914296880563062937443011799707142015647962685071743064933626689093956215640185694480670389308276588763735300266054342210518510655503938301168859241232366499082115661704416570624348536949976639770352500574819937312763118260737810176902179509516452097179033706353327008781474150697920827389789734610788578490948322681570725850626191531687621385181189566240487471233234659263652702257011861689546257355290042178376098634268810888992322511302018026179564203024011436332852715101641086471023435541325718482408540095488497904963", 10).unwrap();
            q = Integer::from_str_radix("3437698264967377277030634390176995339108761179197266750352568811124733920898558792968942018331196975819629005205008386113587643859808501473991043603423893930344847276501635656115934929863928487189206920868313713647633563058406832931496433057597415655426803516291613870104861183026912877156823405411321004129466972852732763205585115072430483394506661540167036456924629674894206529803341592142340709693772564216807751480082209883528322329167890554845107315685301749569796946314108563426847831535198565314370173512451464243428810143307683852607035244704650171763912281349917627223963810996092135117877342840690507979526886087476901475153679297045764348723938821743929206668894364791858405008648706877501022484609478488897352403273816492628362528291051657763913142495017627088665776991501665058363019172356379867940793605694439311020508346728498443971266879607174176083764814467005057407301185462205210470795534594752168817506103", 10).unwrap();
        } else {
            p = Integer::from_str_radix("863", 10).unwrap();
            q = Integer::from_str_radix("859", 10).unwrap();
        }
        let n = p.clone() * q.clone();
        let n_sq = n.clone().square();
        let n_root = n.clone().root(2);
        let mut g = n_sq.clone();
        let mut rand = RandState::new();
        g = g.clone().random_below(&mut rand);
        
        let phi_n: Integer = (p.clone() - 1) * (q.clone() - 1);
        let delta = n.clone() * phi_n.clone();

        // for lambda = g^phi_n, define L = (lambda-1) / n and mu = L^{-1} mod n
        let lambda = g.clone().pow_mod(&phi_n, &n_sq).unwrap();
        let l: Integer = (lambda.clone() - Integer::from(1)) / n.clone();
        let mu = l.clone().invert(&n).unwrap() * phi_n.clone();

        Group {
            n,
            p,
            q,
            n_sq,
            n_root,
            g,
            mu,
            phi_n,
            delta,
            // ... initialize other fields
        }
    }
}

use std::fmt;

impl fmt::Display for Group {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Group:\n")?;
        write!(f, "    n: {}\n", self.n)?;
        write!(f, "    p: {}\n", self.p)?;
        write!(f, "    q: {}\n", self.q)?;
        write!(f, "    n_sq: {}\n", self.n_sq)?;
        write!(f, "    n_root: {}\n", self.n_root)?;
        write!(f, "    g: {}\n", self.g)?;
        write!(f, "    mu: {}\n", self.mu)?;
        write!(f, "    phi_n: {}\n", self.phi_n)?;
        write!(f, "    delta: {}\n", self.delta)?;
        // ... include other fields you want to display

        Ok(())
    }
}

pub fn discrete_logarithm(x: Integer, grp: &Group) -> Integer {
    (x - 1) / grp.n.clone()
}