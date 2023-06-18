// keys.rs

extern crate rug;
use rug::Integer;
use rug::rand::RandState;

use crate::util::{group::Group, matrix::*};
use crate::util::vector::vec_mod;

pub struct IpeSk {
    pub a: Vec<Integer>,
    pub u: Matrix,
    pub u_t: Matrix,
    pub d: Matrix,
    pub d_inv: Matrix,
    pub d_perp: Matrix,
    pub sk1: Vec<Integer>,
    pub sk_enc: Matrix,
}

impl IpeSk {
    pub fn new(dim: usize, b: usize, grp: &Group, rng: &mut RandState<'_>) -> IpeSk {
        // modulus = grp.delta
        // a = random vector of size 2
        // U = random matrix of size dim x 2
        // U_t = transpose of U
        // D = random matrix of size (dim + 2) x (dim + 2 + B)
        // D_inv = right inverse of D
        // D_perp = null space basis of D
        // sk1 = (D_inv_left + D_inv_right * U) * a
        // sk_enc = D_inv_right

        let mut a: Vec<Integer> = Vec::new();
        for _ in 0..2 {
            a.push(grp.delta.clone().random_below(rng));
        }
        // a.push(Integer::from_str_radix("637399315936", 10).unwrap());
        // a.push(Integer::from_str_radix("590061228855", 10).unwrap());
        
    //     let u_parse = "[
    //         [ 417662115945  616437120831 ]
    //         [ 32771204829  27422471670 ]
    //         [ 72208520169  214120036458 ]
    //         [ 720584100637  514501682720 ]
    //       ]";

    //   let d_parse = "[
    //     [ 178153782505  58692816769  545320539919  663932034195  425940656621  201028206110  364737761619  405590618619  9323228293 ]
    //     [ 123343237555  713460155795  450332376329  262524804946  573265876894  246262896207  253288657183  256726918686  531478490427 ]
    //     [ 637591905533  652563464958  131034989664  489643797671  100619175362  246771284412  282711669484  135537507214  86594162909 ]
    //     [ 338768725440  682233745705  2875745563  248754982038  367416689439  344657796082  124806090792  282535386393  454441110691 ]
    //     [ 64997967773  666965197566  181436435360  696822445597  500864376304  542061806807  339033463655  374283287992  699204850384 ]
    //     [ 318272958149  519136877502  528545420885  704131474262  447166133954  447390092006  256033142497  49922644757  634374949304 ]
    //   ]";
    //     let d_inv_parse = "[
    //         [ 326291465227  700194546373  342411885058  352172802613  674343658578  503601581682 ]
    //         [ 181757011142  311142490271  22680357170  402329702169  18013549841  491495038966 ]
    //         [ 679305627469  267733414221  495346811097  220883751166  519548471414  512902545196 ]
    //         [ 234716570367  333365295417  370151620011  99082347247  301175900907  667816846652 ]
    //         [ 166798529756  251510439558  80163735257  171551798428  329487410369  687025062011 ]
    //         [ 698496645543  195798202670  198518718114  47411629245  595946708327  694198851096 ]
    //         [ 432383737699  469279078126  686301861719  50800460781  129040185095  476124616950 ]
    //         [ 596820248085  283406876656  12296429854  463664584186  282316185753  427832870731 ]
    //         [ 444065598142  346411185162  674023788936  433008167936  123287948799  131155170738 ]
    //       ]";
    //     let d_null_parse = "[
    //         [ 138461419988  400323217802  381249536732  700447098017  161996281002  703352966152 ]
    //         [ 478295215476  238386268604  649631871952  516728302453  641249583805  159948011188 ]
    //         [ 684821916476  409358663234  643089502044  489968606582  124666576851  14734255396 ]
    //         [ 654688758631  113071628478  435400700601  482436914680  167795146696  435780740311 ]
    //         [ 422514896576  427402896914  703990092  96606412086  363475573751  468610185568 ]
    //         [ 236057091856  539683574288  620578359252  425047022216  141403826746  547292926708 ]
    //         [ 150987579579  622189929904  417297556765  400587600487  13214063550  432561605811 ]
    //         [ 79872470871  214549392176  331248026421  567994997892  17363158157  169144214887 ]
    //         [ 3625103393  246574003056  329293178291  627989858081  501878198572  294551577405 ]
    //       ]";
        // let u = Matrix::parse_matrix(u_parse);

        let u = Matrix::random(dim, 2, &grp.delta, rng);
        let u_t = u.transpose();
        let (d, d_inv, d_perp) = generate_right_inverse_space(dim + 2, dim + 2 + b, &grp.delta, rng);
        // let d = Matrix::parse_matrix(d_parse);
        // let d_inv = Matrix::parse_matrix(d_inv_parse);
        // let d_perp = Matrix::parse_matrix(d_null_parse);


        // D_inv = (D_inv_left (col=2) || sk_enc (col=dim))

        // sk1 = (D_inv_left + D_inv_right * U) * a
        let mut d_inv_left = Matrix::new(d_inv.rows, 2);
        let mut sk_enc = Matrix::new(d_inv.rows, dim);
        for i in 0..d_inv.rows {
            for j in 0..d_inv.cols {
                if j < 2 {
                    d_inv_left.set(i, j, d_inv.get(i, j));
                } else {
                    sk_enc.set(i, j - 2, d_inv.get(i, j));
                }
            }
        }
        let mut sk1 = (d_inv_left + &(sk_enc.clone() * u.clone())).mul_vec(&a);
        vec_mod(&mut sk1, &grp.delta);
        
        IpeSk {
            a,
            u,
            u_t,
            d,
            d_inv,
            d_perp,
            sk1,
            sk_enc
        }
    }
}

use std::fmt;
impl fmt::Display for IpeSk {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        writeln!(f, "a: {:?}", self.a)?;
        writeln!(f, "u: {:?}", self.u)?;
        writeln!(f, "u_t: {:?}", self.u_t)?;
        writeln!(f, "d: {:?}", self.d)?;
        writeln!(f, "d_inv: {:?}", self.d_inv)?;
        writeln!(f, "d_perp: {:?}", self.d_perp)?;
        writeln!(f, "sk1: {:?}", self.sk1)?;
        writeln!(f, "sk_enc: {:?}", self.sk_enc)?;
        Ok(())
    }
}
