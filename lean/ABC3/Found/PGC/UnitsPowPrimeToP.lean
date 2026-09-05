import ABC3.Found.PGC.UnitsGroupInvariants
import ABC3.Found.PGC.UnitsSplit
import Mathlib.GroupTheory.SpecificGroups.Cyclic
import Mathlib.Data.ZMod.Basic

/-!
# `p ∤ n` のときの単数群・乗法群の `n` 乗指数

[pGC] Proposition 1.2 への**経路 C**(`ResearchPaper/pgc-goal.md`)の計数部分。
経路 C は「`N_n(Γ_K) := #Hom_cont(Γ_K, ℤ/n)` を不変量に取り、`p ∤ n` のとき
`N_n(Γ_K) = n · gcd(n, q−1)` を示す」というものだった。その右辺の由来である
**局所体の乗法群側の指数**を、ここで sorry 無しで確定させる。

## 本ファイルの三つの結論

* **`exists_pow_eq_of_residue_eq_one`**(B1)——主単数(剰余が `1` の単数)は
  `n` 乗根を持つ。`X^n − u` の剰余体での還元は `X^n − 1` で、`1` はその**単根**
  (微分は `n X^{n-1}`、`1` での値 `n` は `p ∤ n` から `𝓀` の中で `≠ 0`)。
  Hensel(`Found/HenselianSplits.lean::exists_root_of_residue_root`)が持ち上げる。
  持ち上げた根の剰余はふたたび `1` なので、**主単数群の中で** `n` 乗写像が全射になる。
* **`index_powRange_units`**(B2)——`[𝒪_K^× : (𝒪_K^×)^n] = gcd(n, q−1)`。
  Teichmüller 持ち上げ(`Found/Teichmuller.lean`)が剰余写像
  `𝒪_K^× ↠ 𝓀^×` の群としての切断を与えるので、`prodKerOfRightInverse` で

  ```
  𝒪_K^×  ≃*  𝓀^× × (1 + 𝔪_K)
  ```

  と分解する。`n` 乗部分群は直積の成分ごとの `n` 乗部分群(`range_powMonoidHom_prod`)
  なので指数は積(`Subgroup.index_prod`)。第二成分は B1 により指数 `1`、
  第一成分は `𝓀^×` が位数 `q−1` の巡回群だから
  `IsCyclic.index_powMonoidHom_range` で `gcd(q−1, n)`。
* **`index_powRange_carrierUnits`**(B3)——`[K^× : (K^×)^n] = n · gcd(n, q−1)`。
  素元を固定した分裂 `K^× ≅ ℤ × 𝒪_K^×`(`Found/PGC/UnitsSplit.lean::unitsSplitEquiv`)
  に B2 と `[ℤ : nℤ] = n`(`index_powRange_multiplicativeInt`)を掛ける。

## ★退化の自己検査——`¬ p ∣ n` を落とすと何が起きるか

`¬ p ∣ n` は装飾ではない。落とすと**主張が偽になる**(自明化ではない)。

* **`n = p` のとき B1 は偽**。`K = ℚ_p`(`p` は奇素数)、`u = 1 + p` を取ると
  `u` の剰余は `1` だが `u` は `ℤ_p^×` の中で `p` 乗数ではない——主単数の `p` 乗は
  `(1+px)^p ∈ 1 + p^2 ℤ_p` に落ちるので、`1 + p ∉ (1 + pℤ_p)^p`。
  Hensel が効かないのは、還元 `X^p − 1 = (X−1)^p` が `1` を**重根**に持ち、
  上の証明で使った `(n : 𝓀) ≠ 0` がちょうど破れるからである。
* **`n = p` のとき B2 も偽**。正しい値は
  `[𝒪_K^× : (𝒪_K^×)^p] = p^{[K:ℚ_p]} · #μ_p(K)`
  (`Found/PGC/UnitsGroupInvariants.lean::index_powRange_smallPrincipalUnits` が
  与える主単数側の `p^{[K:ℚ_p]}` に、`p` 冪根の分 `#μ_p(K)` が掛かる)であって、
  `gcd(p, q−1) = 1` ではない(`q − 1` は `p` と素なので `gcd(p, q−1) = 1`
  ——`PrimeToPTorsion.lean::not_dvd_residueCard_sub_one`)。
  したがって `[K:ℚ_p] ≥ 1` の全ての `K` で右辺と左辺は食い違う。
* **`n = p` のとき B3 も偽**。同様に `p^{[K:ℚ_p]+1} · #μ_p(K) ≠ p · 1`。
* **`n = 0` も除かれている**(`p ∣ 0` なので仮定が排除する)。もし許すと
  `(powMonoidHom 0).range = ⊥` で `𝒪_K^×` は無限だから B2 の左辺は
  `Subgroup.index` の規約で `0`、右辺は `gcd(0, q−1) = q−1 ≠ 0` となり偽。
  (B3 の方は右辺も `0 · (q−1) = 0` になるので、`n = 0` では**偶然に真**
  ——「偽になる」と「自明になる」の両方が同居する珍しい例。)

## 逸脱の記録

原典 [pGC] は Proposition 1.2 の論拠を Serre の局所類体論(相互律)に投げている。
本ファイルはその相互律を経由せず、Hensel と Teichmüller だけで乗法群側の計数を
出している(`ResearchPaper/pgc-goal.md` の「経路 C」で記録済みの逸脱)。
消費する結論は「`q` と `[K:ℚ_p]` が群論的に決まる」だけなので影響は無い。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC ABC3.Found Subgroup Polynomial
open scoped NormedField Valued

variable {p : ℕ} [Fact p.Prime]

/-! ## 準備——`p ∤ n` なら `n` は `𝒪_K` の単数、剰余体でも `≠ 0` -/

/-- `p ∤ n` なら `n` は `𝒪_K` の単数(`‖n‖ = 1`)。 -/
theorem isUnit_natCast_carrierIntegers (K : PAdicLocalField p) {n : ℕ} (hn : ¬ p ∣ n) :
    IsUnit ((n : ℕ) : 𝒪[K.carrier]) := by
  rw [Valued.integer.isUnit_iff_norm_eq_one]
  have hc : ((((n : ℕ) : 𝒪[K.carrier])) : K.carrier) = ((n : ℕ) : K.carrier) := by
    push_cast; ring
  have heq : ‖((n : ℕ) : 𝒪[K.carrier])‖ = ‖((((n : ℕ) : 𝒪[K.carrier])) : K.carrier)‖ := rfl
  rw [heq, hc]
  exact norm_natCast_eq_one_of_not_dvd K hn

/-- したがって剰余体でも `n ≠ 0`——`X^n − 1` が `1` を**単根**に持つことの中身。 -/
theorem natCast_residueField_ne_zero (K : PAdicLocalField p) {n : ℕ} (hn : ¬ p ∣ n) :
    ((n : ℕ) : 𝓀[K.carrier]) ≠ 0 := by
  rw [show ((n : ℕ) : 𝓀[K.carrier])
      = IsLocalRing.residue 𝒪[K.carrier] ((n : ℕ) : 𝒪[K.carrier]) from (map_natCast _ n).symm]
  exact (IsLocalRing.residue_ne_zero_iff_isUnit _).mpr (isUnit_natCast_carrierIntegers K hn)

/-! ## B1——主単数は `n` 乗根を持つ(Hensel) -/

/-- **★★★★主単数群の中で `n` 乗写像は全射**(`p ∤ n`)。

`X^n − u` はモニックで、剰余体での還元 `X^n − 1` は `1` を根に持ち、その微分
`n X^{n-1}` の `1` での値 `n` は `p ∤ n` より `𝓀` の中で `≠ 0`。よって `1` は単根で、
Hensel(`exists_root_of_residue_root`)が `𝒪_K` の根 `a`(剰余 `1`)へ持ち上げる。
`a^n = u` は単元なので `a` も単元。**持ち上げた `v` の剰余がふたたび `1` である**ことが
B2 で効く(主単数群の内部での全射性)。 -/
theorem exists_pow_eq_residue_one_of_residue_eq_one (K : PAdicLocalField p) {n : ℕ}
    (hn : ¬ p ∣ n) (u : (𝒪[K.carrier])ˣ)
    (hu : IsLocalRing.residue 𝒪[K.carrier] (u : 𝒪[K.carrier]) = 1) :
    ∃ v : (𝒪[K.carrier])ˣ, v ^ n = u
      ∧ IsLocalRing.residue 𝒪[K.carrier] (v : 𝒪[K.carrier]) = 1 := by
  haveI := henselianLocalRing_carrierIntegers K
  have hn0 : n ≠ 0 := by rintro rfl; exact hn (dvd_zero p)
  set P : Polynomial 𝒪[K.carrier] := X ^ n - C (u : 𝒪[K.carrier]) with hP
  have hPm : P.Monic := monic_X_pow_sub_C _ hn0
  have hPmap : P.map (IsLocalRing.residue 𝒪[K.carrier]) = X ^ n - 1 := by
    rw [hP]; simp [hu]
  have hroot : Polynomial.eval (1 : 𝓀[K.carrier])
      (P.map (IsLocalRing.residue 𝒪[K.carrier])) = 0 := by
    rw [hPmap]; simp
  have hder : Polynomial.eval (1 : 𝓀[K.carrier])
      (Polynomial.derivative (P.map (IsLocalRing.residue 𝒪[K.carrier]))) ≠ 0 := by
    rw [hPmap]
    simp only [derivative_sub, derivative_X_pow, derivative_one, sub_zero, eval_mul, eval_C,
      eval_pow, eval_X, one_pow, mul_one]
    exact natCast_residueField_ne_zero K hn
  obtain ⟨a, ha, hres⟩ := exists_root_of_residue_root P hPm 1 hroot hder
  have haz : a ^ n = (u : 𝒪[K.carrier]) := by
    have hz : Polynomial.eval a P = 0 := ha
    rw [hP] at hz
    simp only [eval_sub, eval_pow, eval_X, eval_C] at hz
    linear_combination hz
  have hUa : IsUnit a := (isUnit_pow_iff hn0).mp (haz ▸ u.isUnit)
  refine ⟨hUa.unit, Units.ext ?_, ?_⟩
  · rw [Units.val_pow_eq_pow_val, IsUnit.unit_spec, haz]
  · rw [IsUnit.unit_spec]; exact hres

/-- **B1**——剰余が `1` の単数は(`p ∤ n` のとき)`n` 乗根を持つ。 -/
theorem exists_pow_eq_of_residue_eq_one (K : PAdicLocalField p) {n : ℕ} (hn : ¬ p ∣ n)
    (u : (𝒪[K.carrier])ˣ)
    (hu : IsLocalRing.residue 𝒪[K.carrier] (u : 𝒪[K.carrier]) = 1) :
    ∃ v : (𝒪[K.carrier])ˣ, v ^ n = u :=
  let ⟨v, hv, _⟩ := exists_pow_eq_residue_one_of_residue_eq_one K hn u hu
  ⟨v, hv⟩

/-! ## 群論の部品——直積の `n` 乗部分群と `[ℤ : nℤ] = n` -/

/-- 直積の `n` 乗部分群は、成分ごとの `n` 乗部分群の直積。 -/
theorem range_powMonoidHom_prod {G H : Type*} [CommGroup G] [CommGroup H] (n : ℕ) :
    (powMonoidHom n : G × H →* G × H).range
      = ((powMonoidHom n : G →* G).range).prod ((powMonoidHom n : H →* H).range) := by
  ext z
  constructor
  · rintro ⟨⟨x, y⟩, rfl⟩
    exact ⟨⟨x, rfl⟩, ⟨y, rfl⟩⟩
  · rintro ⟨⟨x, hx⟩, ⟨y, hy⟩⟩
    refine ⟨(x, y), ?_⟩
    show ((x, y) : G × H) ^ n = z
    rw [Prod.pow_mk]
    exact Prod.ext hx hy

/-- **`[ℤ : nℤ] = n`**(乗法表記)——`ℤ ↠ ℤ/n` の核が `n` 乗部分群。 -/
theorem index_powRange_multiplicativeInt {n : ℕ} [NeZero n] :
    ((powMonoidHom n : Multiplicative ℤ →* Multiplicative ℤ).range).index = n := by
  set f : Multiplicative ℤ →* Multiplicative (ZMod n) :=
    AddMonoidHom.toMultiplicative (Int.castAddHom (ZMod n)) with hf
  have hfsurj : Function.Surjective f := ZMod.intCast_surjective
  have hker : f.ker = (powMonoidHom n : Multiplicative ℤ →* _).range := by
    ext x
    constructor
    · intro hx
      have h0 : ((Multiplicative.toAdd x : ℤ) : ZMod n) = 0 := hx
      obtain ⟨y, hy⟩ := (ZMod.intCast_zmod_eq_zero_iff_dvd _ n).mp h0
      refine ⟨Multiplicative.ofAdd y, Multiplicative.toAdd.injective ?_⟩
      show n • y = Multiplicative.toAdd x
      rw [hy, nsmul_eq_mul]
    · rintro ⟨y, rfl⟩
      show ((Multiplicative.toAdd (y ^ n) : ℤ) : ZMod n) = 0
      refine (ZMod.intCast_zmod_eq_zero_iff_dvd _ n).mpr ⟨Multiplicative.toAdd y, ?_⟩
      show (n • Multiplicative.toAdd y : ℤ) = (n : ℤ) * Multiplicative.toAdd y
      rw [nsmul_eq_mul]
  rw [← hker]
  show Nat.card (Multiplicative ℤ ⧸ f.ker) = n
  rw [Nat.card_congr (QuotientGroup.quotientKerEquivOfSurjective f hfsurj).toEquiv]
  show Nat.card (ZMod n) = n
  exact Nat.card_zmod n

/-! ## B2——`[𝒪_K^× : (𝒪_K^×)^n] = gcd(n, q−1)` -/

/-- **`𝒪_K^× ≃* 𝓀^× × (1 + 𝔪_K)`**——Teichmüller 持ち上げが剰余写像の切断だから。
第二成分の `(unitsResidueHom K).ker` が主単数群 `1 + 𝔪_K`。 -/
noncomputable def unitsResidueSplit (K : PAdicLocalField p) :
    (𝒪[K.carrier])ˣ ≃* (𝓀[K.carrier])ˣ × (unitsResidueHom K).ker :=
  prodKerOfRightInverse (unitsResidueHom K) (teichmuller 𝒪[K.carrier])
    (residue_teichmuller 𝒪[K.carrier])

/-- **主単数群の中で `n` 乗写像は全射**(B1 の言い換え)。 -/
theorem pow_surjective_principalUnits (K : PAdicLocalField p) {n : ℕ} (hn : ¬ p ∣ n) :
    Function.Surjective
      (powMonoidHom n : (unitsResidueHom K).ker →* (unitsResidueHom K).ker) := by
  rintro ⟨w, hw⟩
  have hres : IsLocalRing.residue 𝒪[K.carrier] (w : 𝒪[K.carrier]) = 1 := by
    have h := congrArg Units.val hw
    simpa using h
  obtain ⟨v, hv, hvres⟩ := exists_pow_eq_residue_one_of_residue_eq_one K hn w hres
  exact ⟨⟨v, Units.ext hvres⟩, Subtype.ext hv⟩

/-- ゆえに主単数群の `n` 乗部分群の指数は `1`——`p ∤ n` の全ての情報がここに入る。 -/
theorem index_powRange_principalUnits (K : PAdicLocalField p) {n : ℕ} (hn : ¬ p ∣ n) :
    ((powMonoidHom n : (unitsResidueHom K).ker →* (unitsResidueHom K).ker).range).index = 1 := by
  rw [Subgroup.index_eq_one]
  exact MonoidHom.range_eq_top_of_surjective _ (pow_surjective_principalUnits K hn)

/-- **★★★★★B2:`[𝒪_K^× : (𝒪_K^×)^n] = gcd(n, q−1)`**(`p ∤ n`)。

`𝒪_K^× ≃* 𝓀^× × (1+𝔪_K)` の第一成分は位数 `q−1` の巡回群で
`IsCyclic.index_powMonoidHom_range` が `gcd(q−1, n)` を与え、第二成分は
Hensel(B1)で指数 `1`。 -/
theorem index_powRange_units (K : PAdicLocalField p) {n : ℕ} (hn : ¬ p ∣ n) :
    ((powMonoidHom n : (𝒪[K.carrier])ˣ →* (𝒪[K.carrier])ˣ).range).index
      = Nat.gcd n (Nat.card 𝓀[K.carrier] - 1) := by
  haveI : DecidableEq 𝓀[K.carrier] := Classical.decEq _
  rw [index_powRange_of_mulEquiv (unitsResidueSplit K) n, range_powMonoidHom_prod,
    Subgroup.index_prod, IsCyclic.index_powMonoidHom_range, index_powRange_principalUnits K hn,
    mul_one, Nat.gcd_comm]
  congr 1
  rw [Nat.card_eq_fintype_card, Fintype.card_units, Nat.card_eq_fintype_card]

/-! ## B3——`[K^× : (K^×)^n] = n · gcd(n, q−1)` -/

/-- **★★★★★B3:`[K^× : (K^×)^n] = n · gcd(n, q−1)`**(`p ∤ n`)。

素元 `ϖ` を一つ取って `K^× ≅ ℤ × 𝒪_K^×`(`unitsSplitEquiv`)。
`ℤ` の側が `n`、`𝒪_K^×` の側が B2 の `gcd(n, q−1)`。
これが経路 C の `N_n(Γ_K) = n · gcd(n, q−1)` の右辺の出どころである。 -/
theorem index_powRange_carrierUnits (K : PAdicLocalField p) {n : ℕ} (hn : ¬ p ∣ n) :
    ((powMonoidHom n : (K.carrier)ˣ →* (K.carrier)ˣ).range).index
      = n * Nat.gcd n (Nat.card 𝓀[K.carrier] - 1) := by
  haveI := valuationRing_isDVR K
  haveI : NeZero n := ⟨by rintro rfl; exact hn (dvd_zero p)⟩
  obtain ⟨ϖ, hϖ⟩ := IsDiscreteValuationRing.exists_irreducible 𝒪[K.carrier]
  rw [← index_powRange_of_mulEquiv (unitsSplitEquiv K hϖ) n, range_powMonoidHom_prod,
    Subgroup.index_prod, index_powRange_multiplicativeInt, index_powRange_units K hn]

end ABC3.Found.PGC
