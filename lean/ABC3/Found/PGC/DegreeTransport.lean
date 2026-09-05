import ABC3.Found.PGC.KummerDuality
import ABC3.Found.PGC.UnitsPowP
import ABC3.Found.PGC.UnitsSplit
import ABC3.Found.PGC.InertiaTransport
import Mathlib.RingTheory.RootsOfUnity.AlgebraicallyClosed

/-!
# 次数の移送——`[K:ℚ_p]` は `Γ_K` の位相群としての同型類だけで決まる

[pGC] Proposition 1.2 への**経路 C**(`ResearchPaper/pgc-goal.md`)の **(C-d)** 本体。
到達点は

  `finrank_eq_of_absGal_equiv :
     ContinuousMulEquiv Γ_K Γ_{K'} → [K:ℚ_p] = [K':ℚ_p]`

である。経路 C は原典の論拠(局所類体論の相互律)を経由しない——その逸脱は
`ResearchPaper/pgc-goal.md` に記録済み。

## 筋

1. **`μ_p ⊆ F` なら `[F^× : (F^×)^p] = p^{[F:ℚ_p]+2}`**
   (`index_powRange_carrierUnits_p_of_mem_mu`)。
   `unitsSplitEquiv`(`F^× ≅ ℤ × 𝒪_F^×`、`Found/PGC/UnitsSplit.lean`)で
   `[F^× : (F^×)^p] = p · [𝒪_F^× : (𝒪_F^×)^p]`。後者は在庫の
   `index_powRange_units_p`(`Found/PGC/UnitsPowP.lean`)で `p^{[F:ℚ_p]} · #μ_p(F)`。
   仮定 `μ_p ⊆ F` から `#μ_p(F) = p` なので合わせて `p^{[F:ℚ_p]+2}`。
2. **Kummer 双対で `Γ_F` の側へ**(`contHomCard_absGal_p_of_mem_mu`)。
   在庫 `card_contHom_eq_index_powRange`(`Found/PGC/KummerDuality.lean`)が
   `#Hom_cont(Γ_F, ℤ/p) = [F^× : (F^×)^p]` を与える。
3. **一般の `K` へ戻す**(`finrank_eq_of_absGal_equiv`)。下記。

## ★H5 の設計——特性的開部分群は使わない(循環の回避)

初期案は `ker(Γ_K → Γ_K^{ab}/(p−1))` という**特性的**開部分群に落とすものだった。
これは循環を呼ぶ:その部分群が**開**であることを言うのに `Hom_cont` の有限性、
すなわち (C-q) の上界が要り、(C-d) が (C-q) に依存してしまう。

本ファイルは **canonical でない開部分群**で済ませる。`α : Γ_K ≃ₜ* Γ_{K'}` を与えられたら、
`K̄` の原始 `p` 乗根 `ζ`・`K̄'` の原始 `p` 乗根 `ζ'` を取って

* `A := (K(ζ)).fixingSubgroup`(開——`IntermediateField.fixingSubgroup_isOpen`)
* `B := α^{-1}((K'(ζ')).fixingSubgroup)`(開——`isOpen_comap_of_continuousMulEquiv`)
* `H := A ⊓ B`、`H' := α(H)`

と置く。`H ≤ A` なので `fixedField H ∋ ζ`、`H' ≤ (K'(ζ')).fixingSubgroup` なので
`fixedField H' ∋ ζ'`。よって 2 を `L_H` と `L_{H'}` に当てられる。
`H.index = H'.index` は `α` が**群同型**であることだけから出る
(`Subgroup.index_map_of_bijective`)——`H` の canonical 性は一度も使わない。

あとは在庫を繋ぐ:

* `absGalFixedFieldCME`(`Found/PGC/AdjoinFieldClosure.lean`、`Γ_{L_H} ≃ₜ* H`)
* `subgroupMapCME`(`Found/PGC/InertiaTransport.lean`、`α` の部分群への制限)
* `finrank_fixedField_eq_index`(`Found/PGC/SubgroupCorrespondenceConstruction.lean`)+
  塔の乗法性 ⟹ `[L_H : ℚ_p] = [Γ_K : H] · [K : ℚ_p]`

`p^{[L_H:ℚ_p]+2} = p^{[L_{H'}:ℚ_p]+2}` から `[L_H:ℚ_p] = [L_{H'}:ℚ_p]`、
すなわち `H.index · [K:ℚ_p] = H.index · [K':ℚ_p]` となり、`H.index > 0` で割れる。

## ★`p = 2` を素通しできる理由

`ζ = −1` が常に体の中にあるので `IsPrimitiveRoot ζ 2` は無条件に取れる。
本ファイルは `X^n − a` の既約性判定(mathlib は奇数 `n` 限定)を**どこでも使わない**
——Kummer 側は `Found/PGC/KummerDuality.lean` が素手の Kummer 指標で組んである。

## ★退化の自己検査——仮定を落とすと何が起きるか

* **`contHomCard_absGal_p_of_mem_mu` の `IsPrimitiveRoot ζ p` を落とすと偽になる**
  (自明化ではない)。`p` を奇素数、`F = ℚ_p` を取る。`[F:ℚ_p] = 1` なので
  式は `p^3` を主張する。しかし実際は
  `[ℚ_p^× : (ℚ_p^×)^p] = p · [ℤ_p^× : (ℤ_p^×)^p] = p · (p^1 · #μ_p(ℤ_p)) = p · p · 1 = p^2`
  ——`μ_p(ℚ_p) = {1}`(`p` 奇のとき `ℤ_p^× ≅ ℤ/(p−1) × ℤ_p` に `p` 捩れは無い)だからで、
  `p^2 ≠ p^3`。仮定はちょうどこの `#μ_p(F)` の因子を `p` に固定するために効いている。
* **`finrank_eq_of_absGal_equiv` の `H` が開であることは必須**。
  `absGalFixedFieldCME` も `fixedFieldLocalField` も開性を引数に取るので、
  落とすと**型が付かない**(`L_H` が `ℚ_p` 上有限次であることが言えない)。
* **`α` が位相同型であることも必須**。群同型だけでは `α(A)` の開性が出ないので、
  上の `H'` が開部分群にならない。ただし `H.index = H'.index` の一手だけは
  群同型で足りる(`Subgroup.index_map_of_bijective`)——そこが本ファイルの設計の要。
* **`p` を一般の `n` に替えても述べられるが、内容が変わる**。`p ∤ n` の側は
  `Found/PGC/UnitsPowPrimeToP.lean` が別の値(`n · gcd(n, q−1)`)を与える。
  経路 C は `p` 側で `[K:ℚ_p]`、`p ∤ n` 側で `q` を読む設計である。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC ABC3.Found Subgroup
open scoped NormedField Valued

variable {p : ℕ} [Fact p.Prime]

/-! ## 一般の群論——直積と `Multiplicative ℤ` の `n` 乗指数 -/

/-- 直積の `n` 乗部分群は成分ごとの `n` 乗部分群の直積。 -/
theorem powRange_prod {G H : Type*} [CommGroup G] [CommGroup H] (n : ℕ) :
    (powMonoidHom n : G × H →* G × H).range
      = ((powMonoidHom n : G →* G).range).prod ((powMonoidHom n : H →* H).range) := by
  ext x
  constructor
  · rintro ⟨y, rfl⟩
    exact ⟨⟨y.1, rfl⟩, ⟨y.2, rfl⟩⟩
  · rintro ⟨⟨a, ha⟩, ⟨b, hb⟩⟩
    exact ⟨(a, b), Prod.ext ha hb⟩

/-- したがって `n` 乗指数は直積で掛け算になる。 -/
theorem index_powRange_prod {G H : Type*} [CommGroup G] [CommGroup H] (n : ℕ) :
    ((powMonoidHom n : G × H →* G × H).range).index
      = ((powMonoidHom n : G →* G).range).index * ((powMonoidHom n : H →* H).range).index := by
  rw [powRange_prod, Subgroup.index_prod]

/-- `[ℤ : nℤ] = n`——`ℤ ↠ ℤ/n` の核がちょうど `n` 乗部分群。 -/
theorem index_powRange_multiplicative_int (n : ℕ) [NeZero n] :
    ((powMonoidHom n : Multiplicative ℤ →* Multiplicative ℤ).range).index = n := by
  set f : Multiplicative ℤ →* Multiplicative (ZMod n) :=
    AddMonoidHom.toMultiplicative (Int.castAddHom (ZMod n)) with hf
  have hfsurj : Function.Surjective f := by
    intro y
    exact ⟨Multiplicative.ofAdd ((Multiplicative.toAdd y).val : ℤ), by
      show ((((Multiplicative.toAdd y).val : ℤ) : ZMod n)) = Multiplicative.toAdd y
      push_cast
      simp⟩
  have hker : f.ker = (powMonoidHom n : Multiplicative ℤ →* Multiplicative ℤ).range := by
    ext x
    constructor
    · intro hx
      have h0 : ((Multiplicative.toAdd x : ℤ) : ZMod n) = 0 := hx
      obtain ⟨k, hk⟩ := (ZMod.intCast_zmod_eq_zero_iff_dvd _ n).mp h0
      refine ⟨Multiplicative.ofAdd k, Multiplicative.toAdd.injective ?_⟩
      show n • k = Multiplicative.toAdd x
      rw [nsmul_eq_mul, ← hk]
    · rintro ⟨y, rfl⟩
      show ((Multiplicative.toAdd ((powMonoidHom n) y) : ℤ) : ZMod n) = 0
      have hy : Multiplicative.toAdd ((powMonoidHom n) y) = n • Multiplicative.toAdd y := rfl
      rw [hy, nsmul_eq_mul]
      push_cast
      simp
  rw [← hker]
  show Nat.card (Multiplicative ℤ ⧸ f.ker) = n
  rw [Nat.card_congr (QuotientGroup.quotientKerEquivOfSurjective f hfsurj).toEquiv]
  show Nat.card (ZMod n) = n
  exact Nat.card_zmod n

/-! ## `μ_p ⊆ F` のときの `p` 捩れ -/

/-- 1 の冪根はノルム 1、したがって `𝒪_F` の中にある。 -/
theorem valuation_le_one_of_isPrimitiveRoot (F : PAdicLocalField p) {ζ : F.carrier}
    (hζ : IsPrimitiveRoot ζ p) : Valued.v ζ ≤ 1 := by
  have hp0 : p ≠ 0 := (Fact.out : p.Prime).ne_zero
  have h1 : ‖ζ‖ ^ p = 1 := by rw [← norm_pow, hζ.pow_eq_one, norm_one]
  have h2 : ‖ζ‖ = 1 := by
    rcases pow_eq_one_iff_cases.mp h1 with h | h | h
    · exact absurd h hp0
    · exact h
    · have hnn := norm_nonneg ζ
      rw [h.1] at hnn
      linarith
  have h : Valued.v ζ = (‖ζ‖₊ : NNReal) := NNReal.eq rfl
  rw [h]
  exact_mod_cast le_of_eq h2

/-- **`μ_p ⊆ F` なら `𝒪_F^×` の `p` 捩れはちょうど `p` 個**。

`ζ` は `𝒪_F` の中にあり(上)、そこで原始 `p` 乗根のままなので
`rootsOfUnity p 𝒪_F = ⟨ζ⟩` は位数 `p` の巡回群。 -/
theorem card_torsion_units_p_of_mem_mu (F : PAdicLocalField p) {ζ : F.carrier}
    (hζ : IsPrimitiveRoot ζ p) :
    Nat.card { u : (𝒪[F.carrier])ˣ // u ^ p = 1 } = p := by
  have hp0 : p ≠ 0 := (Fact.out : p.Prime).ne_zero
  haveI : NeZero p := ⟨hp0⟩
  have hv := valuation_le_one_of_isPrimitiveRoot F hζ
  have hζO : IsPrimitiveRoot (⟨ζ, hv⟩ : 𝒪[F.carrier]) p :=
    IsPrimitiveRoot.of_map_of_injective
      (f := (algebraMap 𝒪[F.carrier] F.carrier : 𝒪[F.carrier] →+* F.carrier)) hζ
      Subtype.val_injective
  have hunit := hζO.isUnit hp0
  have hζu : IsPrimitiveRoot hunit.unit p := hζO.isUnit_unit hp0
  have hzp : Subgroup.zpowers hunit.unit = rootsOfUnity p 𝒪[F.carrier] := hζu.zpowers_eq
  have hcard : Nat.card ↥(rootsOfUnity p 𝒪[F.carrier]) = p := by
    rw [← hzp, Nat.card_zpowers, ← hζu.eq_orderOf]
  exact (Nat.card_congr
    (Equiv.subtypeEquivRight (fun u => (mem_rootsOfUnity p u).symm))).trans hcard

/-! ## (C-d) の計数——`μ_p ⊆ F` のとき `N_p(Γ_F) = p^{[F:ℚ_p]+2}` -/

/-- **★★★★★★★★`μ_p ⊆ F` なら `[F^× : (F^×)^p] = p^{[F:ℚ_p]+2}`**。

`F^× ≅ ℤ × 𝒪_F^×`(`unitsSplitEquiv`)で `p · [𝒪_F^× : (𝒪_F^×)^p]`、
在庫の `index_powRange_units_p` で `p · p^{[F:ℚ_p]} · #μ_p(F)`、
仮定から `#μ_p(F) = p`。 -/
theorem index_powRange_carrierUnits_p_of_mem_mu (F : PAdicLocalField p)
    {ζ : F.carrier} (hζ : IsPrimitiveRoot ζ p) :
    ((powMonoidHom p : (F.carrier)ˣ →* (F.carrier)ˣ).range).index
      = p ^ (Module.finrank ℚ_[p] F.carrier + 2) := by
  haveI : NeZero p := ⟨(Fact.out : p.Prime).ne_zero⟩
  obtain ⟨ϖ, hϖ, ⟨e⟩⟩ := exists_unitsSplitEquiv F
  rw [← index_powRange_of_mulEquiv e p, index_powRange_prod,
    index_powRange_multiplicative_int p, index_powRange_units_p F,
    card_torsion_units_p_of_mem_mu F hζ]
  ring

/-- **★★★★★★★★★★(C-d) 本体——`μ_p ⊆ F` なら `N_p(Γ_F) = p^{[F:ℚ_p]+2}`**。

Kummer 双対(`Found/PGC/KummerDuality.lean::card_contHom_eq_index_powRange`)で
`Γ_F` の連続指標の個数を `[F^× : (F^×)^p]` に移し、上の計数を当てる。

★`IsPrimitiveRoot ζ p` を落とすと**偽になる**(モジュール docstring の自己検査を参照)。 -/
theorem contHomCard_absGal_p_of_mem_mu (F : PAdicLocalField p)
    {ζ : F.carrier} (hζ : IsPrimitiveRoot ζ p) :
    contHomCard F.absGal p = p ^ (Module.finrank ℚ_[p] F.carrier + 2) := by
  haveI : CharZero F.carrier :=
    charZero_of_injective_algebraMap (algebraMap ℚ_[p] F.carrier).injective
  rw [card_contHom_eq_index_powRange F.carrier (Fact.out : p.Prime).ne_zero hζ]
  exact index_powRange_carrierUnits_p_of_mem_mu F hζ

/-! ## `μ_p` を含む開部分群の固定体 -/

/-- 代数閉包には原始 `p` 乗根がある(標数 0)。 -/
theorem exists_isPrimitiveRoot_closure (K : PAdicLocalField p) :
    ∃ ζ : K.closure, IsPrimitiveRoot ζ p := by
  haveI : CharZero K.carrier :=
    charZero_of_injective_algebraMap (algebraMap ℚ_[p] K.carrier).injective
  haveI : NeZero p := ⟨(Fact.out : p.Prime).ne_zero⟩
  haveI : NeZero ((p : ℕ) : K.carrier) := NeZero.charZero
  exact HasEnoughRootsOfUnity.exists_primitiveRoot K.closure p

/-- `K(ζ)` は `K` 上有限次なので、その固定部分群は開。 -/
theorem isOpen_fixingSubgroup_adjoin_singleton (K : PAdicLocalField p) (ζ : K.closure) :
    IsOpen (((IntermediateField.adjoin K.carrier ({ζ} : Set K.closure)).fixingSubgroup :
      Subgroup K.absGal) : Set K.absGal) := by
  have hint : IsIntegral K.carrier ζ := IsAlgebraic.isIntegral (Algebra.IsAlgebraic.isAlgebraic _)
  haveI : FiniteDimensional K.carrier
      (IntermediateField.adjoin K.carrier ({ζ} : Set K.closure)) :=
    IntermediateField.adjoin.finiteDimensional hint
  exact IntermediateField.fixingSubgroup_isOpen _

/-- **`K(ζ)` の固定部分群に含まれる開部分群の固定体は `μ_p` を含む**。

`H ≤ (K(ζ)).fixingSubgroup` なら `H` の元は `ζ` を止めるので `ζ ∈ fixedField H`。
`ζ` は部分体の中でも原始 `p` 乗根のままである(単射準同型で戻す)。 -/
theorem exists_isPrimitiveRoot_fixedFieldLocalField (K : PAdicLocalField p)
    {ζ : K.closure} (hζ : IsPrimitiveRoot ζ p)
    {H : Subgroup K.absGal} (hH : IsOpen (H : Set K.absGal))
    (hle : H ≤ (IntermediateField.adjoin K.carrier ({ζ} : Set K.closure)).fixingSubgroup) :
    ∃ η : (fixedFieldLocalField K H hH).carrier, IsPrimitiveRoot η p := by
  have hmem : ζ ∈ IntermediateField.fixedField H := by
    rw [IntermediateField.mem_fixedField_iff]
    intro f hf
    have h1 := hle hf
    rw [IntermediateField.mem_fixingSubgroup_iff] at h1
    exact h1 ζ (IntermediateField.subset_adjoin K.carrier _ rfl)
  refine ⟨⟨ζ, hmem⟩, ?_⟩
  exact IsPrimitiveRoot.of_map_of_injective
    (f := (IntermediateField.fixedField H).subtype) hζ (fun a b h => Subtype.ext h)

/-- **開部分群の固定体の `ℚ_p` 上の次数**——塔 `ℚ_p ⊆ K ⊆ L_H` の乗法性と
`[L_H : K] = [Γ_K : H]`(`finrank_fixedField_eq_index`)。 -/
theorem finrank_fixedFieldLocalField (K : PAdicLocalField p) (H : Subgroup K.absGal)
    (hH : IsOpen (H : Set K.absGal)) :
    Module.finrank ℚ_[p] (fixedFieldLocalField K H hH).carrier
      = H.index * Module.finrank ℚ_[p] K.carrier := by
  have h : Module.finrank ℚ_[p] K.carrier
        * Module.finrank K.carrier ↥(IntermediateField.fixedField H)
      = Module.finrank ℚ_[p] ↥(IntermediateField.fixedField H) :=
    Module.finrank_mul_finrank ℚ_[p] K.carrier ↥(IntermediateField.fixedField H)
  rw [finrank_fixedField_eq_index K H hH] at h
  show Module.finrank ℚ_[p] ↥(IntermediateField.fixedField H) = _
  rw [← h]
  ring

/-! ## ★(C-d) の到達点——次数の移送 -/

/-- **★★★★★★★★★★★★★★★★★★★★`[K:ℚ_p]` は `Γ_K` の位相群としての
同型類だけで決まる**([pGC] Proposition 1.2 の次数側)。

`α : Γ_K ≃ₜ* Γ_{K'}` から、両辺の `μ_p` を止める開部分群の共通部分 `H`(と `α(H)`)へ
落として `contHomCard_absGal_p_of_mem_mu` を当てる。`H` は canonical に選んでいない
——`H.index = α(H).index` に必要なのは `α` が群同型であることだけなので、それで足りる
(モジュール docstring の「H5 の設計」を参照。特性的部分群を使うと (C-q) と循環する)。 -/
theorem finrank_eq_of_absGal_equiv (K K' : PAdicLocalField p)
    (α : ContinuousMulEquiv K.absGal K'.absGal) :
    Module.finrank ℚ_[p] K.carrier = Module.finrank ℚ_[p] K'.carrier := by
  obtain ⟨ζ, hζ⟩ := exists_isPrimitiveRoot_closure K
  obtain ⟨ζ', hζ'⟩ := exists_isPrimitiveRoot_closure K'
  set A : Subgroup K.absGal :=
    (IntermediateField.adjoin K.carrier ({ζ} : Set K.closure)).fixingSubgroup with hAdef
  set B : Subgroup K.absGal :=
    Subgroup.comap α.toMulEquiv.toMonoidHom
      ((IntermediateField.adjoin K'.carrier ({ζ'} : Set K'.closure)).fixingSubgroup) with hBdef
  have hA : IsOpen (A : Set K.absGal) := isOpen_fixingSubgroup_adjoin_singleton K ζ
  have hB : IsOpen (B : Set K.absGal) :=
    isOpen_comap_of_continuousMulEquiv α _ (isOpen_fixingSubgroup_adjoin_singleton K' ζ')
  have hH : IsOpen ((A ⊓ B : Subgroup K.absGal) : Set K.absGal) := by
    rw [Subgroup.coe_inf]
    exact hA.inter hB
  set H : Subgroup K.absGal := A ⊓ B with hHdef
  set H' : Subgroup K'.absGal := H.map α.toMulEquiv.toMonoidHom with hH'def
  have hH' : IsOpen ((H' : Subgroup K'.absGal) : Set K'.absGal) :=
    isOpen_map_of_continuousMulEquiv α H hH
  obtain ⟨η, hη⟩ := exists_isPrimitiveRoot_fixedFieldLocalField K hζ hH inf_le_left
  have hle' : H' ≤ (IntermediateField.adjoin K'.carrier ({ζ'} : Set K'.closure)).fixingSubgroup := by
    rintro _ ⟨x, hx, rfl⟩
    exact hx.2
  obtain ⟨η', hη'⟩ := exists_isPrimitiveRoot_fixedFieldLocalField K' hζ' hH' hle'
  have cme : ContinuousMulEquiv (fixedFieldLocalField K H hH).absGal
      (fixedFieldLocalField K' H' hH').absGal :=
    ((absGalFixedFieldCME K H hH).trans (subgroupMapCME α H)).trans
      (absGalFixedFieldCME K' H' hH').symm
  have hcard := contHomCard_congr cme p
  rw [contHomCard_absGal_p_of_mem_mu _ hη, contHomCard_absGal_p_of_mem_mu _ hη'] at hcard
  have hdeg : Module.finrank ℚ_[p] (fixedFieldLocalField K H hH).carrier
      = Module.finrank ℚ_[p] (fixedFieldLocalField K' H' hH').carrier := by
    have h2 := Nat.pow_right_injective (Fact.out : p.Prime).two_le hcard
    omega
  rw [finrank_fixedFieldLocalField K H hH, finrank_fixedFieldLocalField K' H' hH'] at hdeg
  have hidx : H'.index = H.index :=
    Subgroup.index_map_of_bijective α.toMulEquiv.bijective H
  rw [hidx] at hdeg
  have hpos : 0 < H.index := by
    rw [← finrank_fixedField_eq_index K H hH]
    haveI := finiteDimensional_fixedField_of_isOpen K H hH
    exact Module.finrank_pos
  exact Nat.eq_of_mul_eq_mul_left hpos hdeg

end ABC3.Found.PGC
