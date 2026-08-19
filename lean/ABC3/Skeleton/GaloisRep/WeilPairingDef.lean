import ABC3.Skeleton.GaloisRep.WeilRoot

/-!
# スケルトン —— **Weil 対 `e_n` とその 4 性質**(`Skeleton`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★層 3–4 —— 対そのもの

`WeilRoot.lean` の `g_P` から

    e_n(P, Q) := τ_Q(g_P) / g_P

で定める。★これが定数(`F` の元)になるのは、`n` 乗が 1 だからである。

## ★★★★消費される 4 性質

(G5) の `det_cyclotomic` が実際に使うのはこの 4 つである:

| 節点 | 内容 |
|---|---|
| `weilPairing_pow_eq_one` | ★★`e_n(P,Q)^n = 1`(`μ_n` に入る) |
| `weilPairing_add_left` | ★★★双線型性 |
| `weilPairing_self` | ★★交代性 `e_n(P,P) = 1` |
| `weilPairing_nondeg` | ★★★★非退化性 |
| `weilPairing_galois` | ★★★★Galois 同変性 |

★★★このうち **Galois 同変性が `det ρ = 円分指標` を出す**——
`σ(e_n(P,Q)) = e_n(σP, σQ) = e_n(P,Q)^{det ρ(σ)}` が `∧²T_l E ≅ ℤ_l(1)` の内容である。
-/

namespace ABC3.Skeleton.GaloisRep

open ABC3.Meta WeierstrassCurve WeierstrassCurve.Affine

variable {F : Type} [Field F] [DecidableEq F]

/-! ## ★★★★★対そのもの -/

/-- ★★★★★**Weil 対** `e_n : E[n] × E[n] → F`。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`WeilRoot.lean` の `g_P` から `e_n(P,Q) = τ_Q(g_P)/g_P` で定める。 -/
noncomputable def weilPairing (W : WeierstrassCurve.Affine F) (n : ℕ) :
    W.Point → W.Point → F :=
  sorry

/-! ## ★★4 性質 -/

/-- ★★**`e_n(P,Q)` は 1 の `n` 乗根である**。 -/
theorem weilPairing_pow_eq_one (W : WeierstrassCurve.Affine F) (n : ℕ) (hn : 1 ≤ n)
    (P Q : W.Point) (hP : n • P = 0) (hQ : n • Q = 0) :
    (weilPairing W n P Q) ^ n = 1 := by
  sorry

/-- ★★★**双線型性**(第 1 変数)。 -/
theorem weilPairing_add_left (W : WeierstrassCurve.Affine F) (n : ℕ) (hn : 1 ≤ n)
    (P₁ P₂ Q : W.Point) (hP₁ : n • P₁ = 0) (hP₂ : n • P₂ = 0) (hQ : n • Q = 0) :
    weilPairing W n (P₁ + P₂) Q = weilPairing W n P₁ Q * weilPairing W n P₂ Q := by
  sorry

/-- ★★**交代性** `e_n(P,P) = 1`。 -/
theorem weilPairing_self (W : WeierstrassCurve.Affine F) (n : ℕ) (hn : 1 ≤ n)
    (P : W.Point) (hP : n • P = 0) :
    weilPairing W n P P = 1 := by
  sorry

/-- ★★★★**非退化性**——`P ≠ 0` なら `e_n(P, ·)` は自明でない。 -/
theorem weilPairing_nondeg (W : WeierstrassCurve.Affine F) (n : ℕ) (hn : 1 ≤ n)
    (P : W.Point) (hP : n • P = 0) (hP0 : P ≠ 0) :
    ∃ Q : W.Point, n • Q = 0 ∧ weilPairing W n P Q ≠ 1 := by
  sorry

/-- ★★★★**Galois 同変性** `σ(e_n(P,Q)) = e_n(σP, σQ)`。

★★これが `det ρ = 円分指標`(`WeilPairing.lean` の `det_galRep_eq_cyclotomic`)を出す。 -/
theorem weilPairing_galois {K L : Type} [Field K] [DecidableEq K] [Field L] [DecidableEq L]
    [Algebra K L] (W : WeierstrassCurve K) (n : ℕ) (σ : L ≃ₐ[K] L)
    (P Q : (W.baseChange L).toAffine.Point) :
    σ (weilPairing (W.baseChange L).toAffine n P Q)
      = weilPairing (W.baseChange L).toAffine n
          (ABC3.Interface.GaloisRep.galPoint W σ P)
          (ABC3.Interface.GaloisRep.galPoint W σ Q) := by
  sorry

/-! ## ★出典の紐付け(`.src`)と、証明が要求するもの(`.needs`) -/

def weilPairing.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対 e_n の定義)",
    sectionId := "genell-thm-3-8" }

def weilPairing.needs : List ProofObligation :=
  [ .otherPaper "GenEll" "Theorem 3.8(Weil 対の構成——g_P の n 乗根の取り出し)" 19,
    .implicitStep
      "★`τ_Q(g_P)/g_P` が **`F` の元である**(定数関数である)こと。n 乗が 1 で、かつ関数体の代数閉包での議論が要る(10-25 ブロック)" 19 ]

def weilPairing_pow_eq_one.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対が 1 の n 乗根に値を取ること)",
    sectionId := "genell-thm-3-8" }

def weilPairing_pow_eq_one.needs : List ProofObligation :=
  [ .otherPaper "GenEll" "Theorem 3.8(Weil 対 e_n の定義)" 19,
    .implicitStep
      "★`g_P(·+Q)^n = f_P([n]·+[n]Q) = f_P([n]·) = g_P(·)^n`——`[n]Q = 0` を使う 1 行だが、引き戻しの合成が可換であることが要る(5-15 ブロック)" 19 ]

def weilPairing_add_left.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の双線型性)",
    sectionId := "genell-thm-3-8" }

def weilPairing_add_left.needs : List ProofObligation :=
  [ .citation "[Silverman]" "The Arithmetic of Elliptic Curves, III.8.1(a)(双線型性)"
      (.absent "mathlib に Weil 対は 0 件(2026-08-20、`WeilPairing|weil_pairing` で全文検索して 0 件)") 19,
    .otherPaper "GenEll" "Theorem 3.8(Weil 対 e_n の定義)" 19,
    .implicitStep
      "★★第 1 変数の加法性は `f_{P₁+P₂} = f_{P₁}·f_{P₂}·(主因子)` に帰着する。主因子の因子が 0 であることが要る(15-40 ブロック)" 19,
    .implicitStep
      "★第 2 変数の加法性は `τ_{Q₁+Q₂} = τ_{Q₁}∘τ_{Q₂}`(平行移動が群作用であること)から出る(5-10 ブロック)" 19 ]

def weilPairing_self.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の交代性)",
    sectionId := "genell-thm-3-8" }

def weilPairing_self.needs : List ProofObligation :=
  [ .citation "[Silverman]" "The Arithmetic of Elliptic Curves, III.8.1(b)(交代性)"
      (.absent "mathlib に Weil 対は 0 件(2026-08-20、同上の検索)") 19,
    .otherPaper "GenEll" "Theorem 3.8(Weil 対 e_n の定義)" 19,
    .implicitStep
      "★★`∏_{j} τ_{jP}(g_P)` の telescoping。`E[n]` 上の積を取る計算である(10-25 ブロック)" 19 ]

def weilPairing_nondeg.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の非退化性)",
    sectionId := "genell-thm-3-8" }

def weilPairing_nondeg.needs : List ProofObligation :=
  [ .citation "[Silverman]" "The Arithmetic of Elliptic Curves, III.8.1(c)(非退化性)"
      (.absent "mathlib に Weil 対は 0 件(2026-08-20、同上の検索)") 19,
    .otherPaper "GenEll" "Theorem 3.8(Weil 対 e_n の定義)" 19,
    .implicitStep
      "★★★`e_n(P,·) ≡ 1` なら `g_P` が `[n]^*` の像に入り、`f_P` が `n` 乗になる——そこから `P = 0` を出す。Kummer 理論型の議論(20-50 ブロック)" 19,
    .implicitStep
      "★`E[n] ≃ (ℤ/n)²`(第 65-72 ブロック、**Found に済**)がここで消費される(0 ブロック)" 19 ]

def weilPairing_galois.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の Galois 同変性)",
    sectionId := "genell-thm-3-8" }

def weilPairing_galois.needs : List ProofObligation :=
  [ .citation "[Silverman]" "The Arithmetic of Elliptic Curves, III.8.1(d)(Galois 同変性)"
      (.absent "mathlib に Weil 対は 0 件(2026-08-20、同上の検索)") 19,
    .otherPaper "GenEll" "Theorem 3.8(Weil 対 e_n の定義)" 19,
    .implicitStep
      "★★構成が関手的なら自動だが、明示式(`τ_Q(g_P)/g_P`)からだと `σ(f_P) = f_{σP}` と `σ(g_P) = g_{σP}` を別途示すことになる(10-25 ブロック)" 19,
    .implicitStep
      "★`galPoint`(σ の点への作用)は **Interface に定義済**であり posit ではない(0 ブロック)" 19 ]

end ABC3.Skeleton.GaloisRep
