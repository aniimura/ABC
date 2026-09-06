import ABC3.Found.PGC.UnramifiedZhat
import ABC3.Found.PGC.TorsionCyclotome
import ABC3.Found.PGC.CyclotomicRecovery

/-!
# 経路 Λ9 —— `tors_{p^n}(𝒪_L^× × Ẑ) ≅ μ_{p^n}(L)`(作用つき、`sorry` 無し)

[pGC] Proposition 1.1 の壁 `TorsionCyclotomeIsCyclotomic`
(`Found/PGC/CyclotomeTransport.lean`)は Λ7f ＋ Λ8 ＋ Λ9 の 3 つに分かれる。
本ファイルはその **Λ9**——局所類体論が与える群

```
Gal(L_π · L^ur / L) ≅ 𝒪_L^× × Ẑ
```

の **`p^n` 捩れが `μ_{p^n}(L)` である**ことを、**作用まで込めて**示す。

## 筋(3 段)

1. **`Ẑ` は捩れ自由**(`zhat_eq_one_of_pow_eq_one`)。
   `Ẑ = lim_H (ℤ/H)` の元 `x` が `x^m = 1`(`m ≠ 0`)を満たすとする。
   有限指数正規部分群 `H ≤ ℤ` を取り、`N := [ℤ:H]`、`H' := ker(ℤ → ℤ/mN)` とすると
   `H' ≤ H`(`Subgroup.pow_index_mem` から `N ∈ H`)。
   `x_{H'} = mk a` と書くと `a^m ∈ H'`、すなわち `mN ∣ m·a`、よって `N ∣ a`、
   よって `a ∈ H`。極限の整合性から `x_H = mk a = 1`。`H` は任意だったので `x = 1`。
   ★在庫調査(2026-09-06):**`Ẑ` の捩れ自由性は mathlib にも ABC3 にも無かった**
   (`.cache/mathlib-index.txt` の `ProfiniteCompletion` 系は普遍性と稠密性だけ)。
   ここで作った。同時に **`CommGroup Ẑ`** も無かったので作った(`zhatCommGroup`)。
2. **`𝒪_L^×` の `m` 捩れ = `μ_m(L)`**(`torsionUnitsZHatEquiv`)。
   ★ここは**単数群の構造定理 `𝒪^× ≅ μ × ℤ_p^d` を経由しない**。
   `x^m = 1`(`m ≠ 0`)なら `‖x‖^m = 1` かつ `‖x‖ ≥ 0` なので `‖x‖ = 1`、
   したがって `x` も `x⁻¹` も `𝒪_L` に入る——これだけで
   `𝒪_L^× → L^×` は `m` 乗根の上への全単射になる。
   構造定理(`Found/PGC/UnitsSplit.lean` / `PrincipalUnitsRank.lean` /
   `UnitsPowP.lean`)は**使わなかった**(見積より安く済んだ理由)。
3. **作用**(`torsionUnitsZHatEquiv_pow` / `torsionUnitsZHatEquiv_galois` /
   `ringEquiv_pow_eq_cyclotomicCharacter`)。

## ★大きさだけでなく作用が乗っていること

段取りの測定が名指しした危険——「捩れの**大きさ**は出るが**作用**が乗らず
Λ8 が接続できない」——を避けるため、本ファイルは同型を作るだけで終わらせない。

* `torsionUnitsZHatEquiv_pow`:`𝒪_L^× × Ẑ` の**任意の**群自己同型 `Φ` が
  捩れの上で `u ↦ u^c` を誘導するなら、同型 `e` を通した像も `ζ ↦ ζ^c`。
  ★`Φ` が 2 つの因子を混ぜてもよい——捩れ元の `Ẑ` 成分が `1` である
  (`torsion_snd_eq_one`、1 の帰結)ので、混ぜても第 1 成分しか見えない。
* `torsionUnitsZHatEquiv_galois`:同じことを体の自己同型 `σ` の形で
  (`Φ` が `σ` を覆うなら `e` は `σ`-同変)。
* `ringEquiv_pow_eq_cyclotomicCharacter`:`K ⊆ F̄` の中で `σ` が `g ∈ Γ_F` の
  制限なら、`μ_{p^n}(K)` の上で `σ ζ = ζ^{χ_{F,n}(g)}`。
  ★これが「作用は円分指標倍」の中身であり、
  `exists_torsionUnitsZHat_equiv_cyclotomic` が 3 つを 1 本に束ねる。
  結論の指数は `TorsionCyclotomeIsCyclotomic` と**同じ式**
  `(PadicInt.toZModPow n (cyclotomicCharacter F.closure p g.toRingEquiv)).val` である。

## ★設計上の注意(守ったこと)

* **`Abelianization` を使っていない**。`Ẑ` は `ProfiniteCompletion` の極限記述の
  中でしか商を作らない。
* **`inertia` を経由していない**——`Interface` の `SubgroupCorrespondence` /
  `ResidueCardinality` は本ファイルの主張にも証明にも現れない。
* **自由なパラメータを結論に出していない**——到達点
  `exists_torsionUnitsZHat_equiv_rootsOfUnity` /
  `exists_torsionUnitsZHat_equiv_cyclotomic` の型は `K`(と `F`・`ι`)にしか
  依存せず、同型 `e` は `∃` の内側にある。
* **`sorry` 本体の `def` を作っていない**(D21)。本ファイルに `sorry` は無い。

## 逸脱(記録)

* **`Ẑ` の与え方**。`Ẑ` は `Found/PGC/UnramifiedZhat.lean` に従い
  `ProfiniteGrp.ProfiniteCompletion.completion (GrpCat.of (Multiplicative ℤ))`
  として扱う(mathlib の `ProfiniteCompletion` が乗法群の圏 `GrpCat` の上にあるため)。
  古典的な `Ẑ = lim ℤ/n` との一致は主張していない——本ファイルが使うのは
  「有限指数正規部分群での極限」という定義そのものだけである。
* **`μ_{p^n}(L)` の与え方**。`rootsOfUnity (p^n) L.carrier`(`Subgroup (L^×)`)を使う。
  `Found/PGC/DegreeTransport.lean` が既に `rootsOfUnity p 𝒪[F.carrier]` を
  使っており、そちらとの往復が本ファイルの `unitsToField` である。
* **作用の入力の形**。原典(局所類体論)は「相互律が Galois 同変である」と
  一言で畳むが、Λ8 がまだ無いので、本ファイルは**同変性を仮説として受け取り
  結論へ運ぶ**形(`hΦ`)で書いた。Λ8 が `Φ` と `hΦ` を供給すれば
  そのまま `TorsionCyclotomeIsCyclotomic` の指数の式が出る。
  ★これは前提の追加ではなく、**接続点の明示**である。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open CategoryTheory ProfiniteGrp ProfiniteGrp.ProfiniteCompletion
open scoped NormedField Valued

/-! ## 1. `Ẑ` の中の段 `kℤ` -/

/-- `Multiplicative ℤ →* Multiplicative (ZMod k)`。核は `kℤ`。 -/
def intModHom (k : ℕ) : Multiplicative ℤ →* Multiplicative (ZMod k) :=
  AddMonoidHom.toMultiplicative (Int.castAddHom (ZMod k))

theorem intModHom_apply (k : ℕ) (x : Multiplicative ℤ) :
    intModHom k x = Multiplicative.ofAdd ((Multiplicative.toAdd x : ℤ) : ZMod k) := rfl

theorem mem_ker_intModHom {k : ℕ} {x : Multiplicative ℤ} :
    x ∈ (intModHom k).ker ↔ (k : ℤ) ∣ Multiplicative.toAdd x := by
  rw [MonoidHom.mem_ker, intModHom_apply, ← ZMod.intCast_zmod_eq_zero_iff_dvd]
  exact ⟨fun h => by simpa using congrArg Multiplicative.toAdd h,
    fun h => by rw [h]; rfl⟩

instance finite_quotient_ker_intModHom (k : ℕ) [NeZero k] :
    Finite (Multiplicative ℤ ⧸ (intModHom k).ker) :=
  Finite.of_equiv _ (QuotientGroup.quotientKerEquivRange (intModHom k)).toEquiv.symm

/-- `kℤ ⊆ ℤ` を有限指数正規部分群として。`Ẑ` の極限の添字に使う。

退化の自己検査:`NeZero k` を落とすと**偽**——`k = 0` では核が `⊥` になり、
`ℤ` の中で無限指数なので `FiniteIndexNormalSubgroup` を作れない。 -/
def levelSubgroup (k : ℕ) [NeZero k] : FiniteIndexNormalSubgroup (Multiplicative ℤ) where
  toSubgroup := (intModHom k).ker
  isFiniteIndex' := Subgroup.finiteIndex_of_finite_quotient

/-! ## 2. `Ẑ` は可換で捩れ自由 -/

/-- 副有限極限の冪は成分ごと。 -/
theorem limit_pow_val {J : Type} [SmallCategory J] (F : J ⥤ ProfiniteGrp.{0})
    (x : (limit F).toProfinite.toTop) (j : J) (k : ℕ) :
    (x ^ k).val j = (x.val j) ^ k := by
  induction k with
  | zero => simp
  | succ k ih => rw [pow_succ, ProfiniteGrp.limit_mul_val, ih, pow_succ]

theorem zhat_pow_val (x : ZHat) (H : FiniteIndexNormalSubgroup (Multiplicative ℤ)) (k : ℕ) :
    (x ^ k).val H = (x.val H) ^ k :=
  limit_pow_val _ x H k

theorem zhat_val_one (H : FiniteIndexNormalSubgroup (Multiplicative ℤ)) :
    ((1 : ZHat)).val H = 1 :=
  ProfiniteGrp.limit_one_val _ H

theorem zhat_mul_comm (x y : ZHat) : x * y = y * x := by
  refine ProfiniteGrp.limit_ext _ (x * y) (y * x) (fun H => ?_)
  have h1 : ((x * y).val H) = (x.val H) * (y.val H) := ProfiniteGrp.limit_mul_val _ x y H
  have h2 : ((y * x).val H) = (y.val H) * (x.val H) := ProfiniteGrp.limit_mul_val _ y x H
  rw [h1, h2]
  exact @mul_comm (Multiplicative ℤ ⧸ H.toSubgroup) _ (x.val H) (y.val H)

/-- `Ẑ` は可換群。★mathlib は `ProfiniteCompletion` に可換性を付けていない
(2026-09-06 在庫調査)ので、ここで付ける。 -/
noncomputable instance zhatCommGroup : CommGroup ZHat :=
  { (inferInstance : Group ZHat) with mul_comm := zhat_mul_comm }

/-- **★★★★★★★★★★`Ẑ` は捩れ自由**。

`x^m = 1`(`m ≠ 0`)なら `x = 1`。

段 `H`(有限指数正規、`N := [ℤ:H]`)ごとに、より深い段
`H' := ker(ℤ → ℤ/mN) = mNℤ` を取る。`N ∈ H`(`Subgroup.pow_index_mem`)なので
`H' ≤ H`。`x_{H'} = mk a` と書くと `a^m ∈ H'`、すなわち `mN ∣ m·a`、
`m ≠ 0` を約して `N ∣ a`、よって `a ∈ H`。極限の整合性から `x_H = mk a = 1`。

退化の自己検査:`m ≠ 0` は落とせない——`Ẑ` は自明でない
(`exists_zhat_ne_one`)ので、`m = 0` では結論が偽になる。 -/
theorem zhat_eq_one_of_pow_eq_one {m : ℕ} (hm : m ≠ 0) {x : ZHat} (hx : x ^ m = 1) : x = 1 := by
  refine ProfiniteGrp.limit_ext _ x 1 (fun H => ?_)
  rw [ProfiniteGrp.limit_one_val]
  haveI := H.isNormal'
  have hN : H.toSubgroup.index ≠ 0 := Subgroup.finiteIndex_iff.mp H.isFiniteIndex'
  haveI : NeZero (m * H.toSubgroup.index) := ⟨Nat.mul_ne_zero hm hN⟩
  have hgen : (Multiplicative.ofAdd ((H.toSubgroup.index : ℕ) : ℤ)) ∈ H.toSubgroup := by
    have h := Subgroup.pow_index_mem H.toSubgroup (Multiplicative.ofAdd (1 : ℤ))
    have key : (Multiplicative.ofAdd ((H.toSubgroup.index : ℕ) : ℤ))
        = Multiplicative.ofAdd (1 : ℤ) ^ H.toSubgroup.index := by
      rw [← ofAdd_nsmul, nsmul_eq_mul, mul_one]
    rw [key]; exact h
  have hle : levelSubgroup (m * H.toSubgroup.index) ≤ H := by
    intro y hy
    obtain ⟨j, hj⟩ := mem_ker_intModHom.mp hy
    have hy' : y = (Multiplicative.ofAdd ((H.toSubgroup.index : ℕ) : ℤ)) ^ ((m : ℤ) * j) := by
      apply Multiplicative.toAdd.injective
      rw [toAdd_zpow, hj]
      simp only [toAdd_ofAdd, smul_eq_mul]
      push_cast
      ring
    rw [hy']
    exact zpow_mem hgen _
  obtain ⟨a, ha⟩ : ∃ a : Multiplicative ℤ,
      (QuotientGroup.mk a :
          Multiplicative ℤ ⧸ (levelSubgroup (m * H.toSubgroup.index)).toSubgroup)
        = ((x.val (levelSubgroup (m * H.toSubgroup.index)) : _) :
            Multiplicative ℤ ⧸ (levelSubgroup (m * H.toSubgroup.index)).toSubgroup) :=
    QuotientGroup.mk_surjective _
  have hpow : ((x.val (levelSubgroup (m * H.toSubgroup.index)) : _) :
      Multiplicative ℤ ⧸ (levelSubgroup (m * H.toSubgroup.index)).toSubgroup) ^ m = 1 := by
    have h1 : (x.val (levelSubgroup (m * H.toSubgroup.index))) ^ m = 1 := by
      rw [← zhat_pow_val, hx, zhat_val_one]
    exact h1
  have hmem : a ^ m ∈ (levelSubgroup (m * H.toSubgroup.index)).toSubgroup := by
    refine (QuotientGroup.eq_one_iff _).mp ?_
    have hrfl : (QuotientGroup.mk (a ^ m) :
          Multiplicative ℤ ⧸ (levelSubgroup (m * H.toSubgroup.index)).toSubgroup)
        = (QuotientGroup.mk a :
          Multiplicative ℤ ⧸ (levelSubgroup (m * H.toSubgroup.index)).toSubgroup) ^ m := rfl
    rw [hrfl, ha]; exact hpow
  have haH : a ∈ H.toSubgroup := by
    obtain ⟨j, hj⟩ := mem_ker_intModHom.mp hmem
    rw [toAdd_pow, nsmul_eq_mul] at hj
    push_cast at hj
    have hja : Multiplicative.toAdd a = (H.toSubgroup.index : ℤ) * j := by
      have hm0 : (m : ℤ) ≠ 0 := Int.natCast_ne_zero.mpr hm
      apply mul_left_cancel₀ hm0
      rw [hj]; ring
    have hfa : a = (Multiplicative.ofAdd ((H.toSubgroup.index : ℕ) : ℤ)) ^ j := by
      apply Multiplicative.toAdd.injective
      rw [toAdd_zpow, hja]
      simp [mul_comm]
    rw [hfa]
    exact zpow_mem hgen _
  have hcompat := x.2 hle.hom
  show ((x.val H : _) : Multiplicative ℤ ⧸ H.toSubgroup) = 1
  have hval : ((x.val H : _) : Multiplicative ℤ ⧸ H.toSubgroup) = QuotientGroup.mk a := by
    rw [← hcompat, ← ha]; rfl
  rw [hval]
  exact (QuotientGroup.eq_one_iff _).mpr haH

theorem zhat_ne_one_eta : (ProfiniteCompletion.etaFn (GrpCat.of (Multiplicative ℤ))
    (Multiplicative.ofAdd (1 : ℤ)) : ZHat) ≠ 1 := by
  intro hcon
  haveI : NeZero 2 := ⟨two_ne_zero⟩
  have h1 : ((ProfiniteCompletion.etaFn (GrpCat.of (Multiplicative ℤ))
      (Multiplicative.ofAdd (1 : ℤ)) : ZHat)).val (levelSubgroup 2) = 1 := by
    rw [hcon]; exact zhat_val_one _
  have h2 : (QuotientGroup.mk (Multiplicative.ofAdd (1 : ℤ)) :
      Multiplicative ℤ ⧸ (levelSubgroup 2).toSubgroup) = 1 := h1
  have h3 := (QuotientGroup.eq_one_iff _).mp h2
  have h4 := mem_ker_intModHom.mp h3
  norm_num at h4

/-- `Ẑ` は自明でない —— `zhat_eq_one_of_pow_eq_one` の `m ≠ 0` は落とせない。 -/
theorem exists_zhat_ne_one : ∃ x : ZHat, x ≠ 1 :=
  ⟨_, zhat_ne_one_eta⟩

/-! ## 3. `m` 捩れ部分群 -/

/-- **`A` の `m` 捩れ部分群**。`Found/PGC/TorsionCyclotome.lean` の `cyclotome` と
同じ `powMonoidHom` の核だが、こちらは位相アーベル化を挟まない。 -/
def powTorsion (A : Type*) [CommGroup A] (m : ℕ) : Subgroup A :=
  (powMonoidHom m : A →* A).ker

@[simp] theorem mem_powTorsion {A : Type*} [CommGroup A] {m : ℕ} {x : A} :
    x ∈ powTorsion A m ↔ x ^ m = 1 := Iff.rfl

theorem powTorsion_fst {A B : Type*} [CommGroup A] [CommGroup B] {m : ℕ} {x : A × B}
    (h : x ∈ powTorsion (A × B) m) : x.1 ^ m = 1 :=
  congrArg Prod.fst h

theorem powTorsion_snd {A B : Type*} [CommGroup A] [CommGroup B] {m : ℕ} {x : A × B}
    (h : x ∈ powTorsion (A × B) m) : x.2 ^ m = 1 :=
  congrArg Prod.snd h

theorem mem_powTorsion_prod {A B : Type*} [CommGroup A] [CommGroup B] {m : ℕ} {a : A} {b : B}
    (ha : a ^ m = 1) (hb : b ^ m = 1) : (a, b) ∈ powTorsion (A × B) m :=
  Prod.ext ha hb

theorem map_powTorsion_of_mulEquiv {A B : Type*} [CommGroup A] [CommGroup B] (e : A ≃* B)
    (m : ℕ) : (powTorsion A m).map e.toMonoidHom = powTorsion B m :=
  map_powKer_of_mulEquiv e m

/-- 群同型が誘導する `m` 捩れの同型。 -/
noncomputable def powTorsionCongr {A B : Type*} [CommGroup A] [CommGroup B] (e : A ≃* B) (m : ℕ) :
    ↥(powTorsion A m) ≃* ↥(powTorsion B m) :=
  (e.subgroupMap (powTorsion A m)).trans (MulEquiv.subgroupCongr (map_powTorsion_of_mulEquiv e m))

@[simp] theorem powTorsionCongr_coe {A B : Type*} [CommGroup A] [CommGroup B] (e : A ≃* B)
    (m : ℕ) (x : ↥(powTorsion A m)) :
    ((powTorsionCongr e m x : ↥(powTorsion B m)) : B) = e (x : A) := rfl

/-! ## 4. `𝒪_K^×` の `m` 捩れ = `μ_m(K)`

★構造定理 `𝒪^× ≅ μ × ℤ_p^d` は使わない。`m` 乗根のノルムが `1` であること
(それ自身と逆元がともに `𝒪_K` に入ること)だけで足りる。 -/

variable {p : ℕ} [Fact p.Prime]

/-- `𝒪_K^× → K^×`。 -/
noncomputable def unitsToField (K : PAdicLocalField p) : (𝒪[K.carrier])ˣ →* (K.carrier)ˣ :=
  Units.map (algebraMap 𝒪[K.carrier] K.carrier : 𝒪[K.carrier] →+* K.carrier).toMonoidHom

@[simp] theorem unitsToField_coe (K : PAdicLocalField p) (u : (𝒪[K.carrier])ˣ) :
    ((unitsToField K u : (K.carrier)ˣ) : K.carrier) = ((u : 𝒪[K.carrier]) : K.carrier) := rfl

theorem unitsToField_injective (K : PAdicLocalField p) : Function.Injective (unitsToField K) := by
  intro a b hab
  exact Units.ext (Subtype.ext (congrArg (fun z : (K.carrier)ˣ => (z : K.carrier)) hab))

/-- `m` 乗根(`m ≠ 0`)はノルム `1`。 -/
theorem carrier_norm_eq_one_of_pow_eq_one (K : PAdicLocalField p) {m : ℕ} (hm : m ≠ 0)
    {x : K.carrier} (hx : x ^ m = 1) : ‖x‖ = 1 :=
  (pow_eq_one_iff_of_nonneg (norm_nonneg x) hm).mp (by rw [← norm_pow, hx, norm_one])

/-- **`m` 乗根は `𝒪_K` に入る**。

退化の自己検査:`m ≠ 0` は落とせない——`m = 0` では `x^0 = 1` が恒真なので、
主張は「`K = 𝒪_K`」になり偽。 -/
theorem mem_carrierIntegers_of_pow_eq_one (K : PAdicLocalField p) {m : ℕ} (hm : m ≠ 0)
    {x : K.carrier} (hx : x ^ m = 1) : x ∈ 𝒪[K.carrier] := by
  show Valued.v x ≤ 1
  have hv : Valued.v x = (‖x‖₊ : NNReal) := NNReal.eq rfl
  rw [hv]
  exact_mod_cast le_of_eq (carrier_norm_eq_one_of_pow_eq_one K hm hx)

/-- `K^×` の `m` 乗根を `𝒪_K^×` の元として読む。逆元も `m` 乗根なので単数になる。 -/
noncomputable def unitOfRootOfUnity (K : PAdicLocalField p) {m : ℕ} (hm : m ≠ 0)
    {y : (K.carrier)ˣ} (hy : y ^ m = 1) : (𝒪[K.carrier])ˣ where
  val := ⟨(y : K.carrier), mem_carrierIntegers_of_pow_eq_one K hm
    (by rw [← Units.val_pow_eq_pow_val, hy, Units.val_one])⟩
  inv := ⟨((y⁻¹ : (K.carrier)ˣ) : K.carrier), mem_carrierIntegers_of_pow_eq_one K hm
    (by rw [← Units.val_pow_eq_pow_val, inv_pow, hy, inv_one, Units.val_one])⟩
  val_inv := Subtype.ext (by simp)
  inv_val := Subtype.ext (by simp)

@[simp] theorem unitsToField_unitOfRootOfUnity (K : PAdicLocalField p) {m : ℕ} (hm : m ≠ 0)
    {y : (K.carrier)ˣ} (hy : y ^ m = 1) : unitsToField K (unitOfRootOfUnity K hm hy) = y :=
  Units.ext rfl

theorem coe_pow_eq_one_of_units_pow_eq_one (K : PAdicLocalField p) {m : ℕ}
    {u : (𝒪[K.carrier])ˣ} (hu : u ^ m = 1) : (((u : 𝒪[K.carrier]) : K.carrier)) ^ m = 1 := by
  have h := congrArg (fun w : (𝒪[K.carrier])ˣ => ((w : 𝒪[K.carrier]) : K.carrier)) hu
  simpa using h

theorem units_eq_pow_of_coe (K : PAdicLocalField p) {u v : (𝒪[K.carrier])ˣ}
    (h : ((unitsToField K u : (K.carrier)ˣ) : K.carrier)
      = ((unitsToField K v : (K.carrier)ˣ) : K.carrier)) : u = v :=
  unitsToField_injective K (Units.ext h)

/-! ## 5. ★到達点(その 1)—— `tors_m(𝒪_K^× × Ẑ) ≅ μ_m(K)` -/

/-- **★★★★★★★★★★★★★★★(Λ9)`tors_m(𝒪_K^× × Ẑ) ≅ μ_m(K)`**。

`Ẑ` が捩れ自由(`zhat_eq_one_of_pow_eq_one`)なので捩れは第 1 因子だけに乗り、
`𝒪_K^×` の `m` 捩れは `K^×` の `m` 捩れ、すなわち `μ_m(K)` に一致する。

退化の自己検査。

* `m ≠ 0` は落とせない——`m = 0` では両辺とも `⊤` になり、主張は
  `𝒪_K^× × Ẑ ≅ K^×` になる(`𝒪_K^×` は素元を含まないので偽)。
* `Ẑ` を捩れのある副有限群(たとえば `ℤ/m`)に替えると**偽**——
  左辺が `m` 倍だけ大きくなる。`zhat_eq_one_of_pow_eq_one` がそこを守っている。 -/
noncomputable def torsionUnitsZHatEquiv (K : PAdicLocalField p) {m : ℕ} (hm : m ≠ 0) :
    ↥(powTorsion ((𝒪[K.carrier])ˣ × ZHat) m) ≃* ↥(rootsOfUnity m K.carrier) where
  toFun x := ⟨unitsToField K (x : (𝒪[K.carrier])ˣ × ZHat).1, by
    rw [mem_rootsOfUnity, ← map_pow, powTorsion_fst x.2, map_one]⟩
  invFun y := ⟨(unitOfRootOfUnity K hm ((mem_rootsOfUnity m _).mp y.2), 1),
    mem_powTorsion_prod (unitsToField_injective K
      (by rw [map_pow, unitsToField_unitOfRootOfUnity, (mem_rootsOfUnity m _).mp y.2, map_one]))
      (one_pow m)⟩
  left_inv x := by
    refine Subtype.ext (Prod.ext ?_ ?_)
    · exact unitsToField_injective K (by rw [unitsToField_unitOfRootOfUnity])
    · exact (zhat_eq_one_of_pow_eq_one hm (powTorsion_snd x.2)).symm
  right_inv y := Subtype.ext (by simp [unitsToField_unitOfRootOfUnity])
  map_mul' a b := Subtype.ext (map_mul (unitsToField K) _ _)

@[simp] theorem torsionUnitsZHatEquiv_coe (K : PAdicLocalField p) {m : ℕ} (hm : m ≠ 0)
    (x : ↥(powTorsion ((𝒪[K.carrier])ˣ × ZHat) m)) :
    ((torsionUnitsZHatEquiv K hm x : ↥(rootsOfUnity m K.carrier)) : (K.carrier)ˣ)
      = unitsToField K (x : (𝒪[K.carrier])ˣ × ZHat).1 := rfl

/-- 捩れ元の `Ẑ` 成分は自明。★作用の同変性が「因子を混ぜる `Φ`」でも通る理由。 -/
theorem torsion_snd_eq_one (K : PAdicLocalField p) {m : ℕ} (hm : m ≠ 0)
    (x : ↥(powTorsion ((𝒪[K.carrier])ˣ × ZHat) m)) :
    (x : (𝒪[K.carrier])ˣ × ZHat).2 = 1 :=
  zhat_eq_one_of_pow_eq_one hm (powTorsion_snd x.2)

theorem torsion_eq_mk (K : PAdicLocalField p) {m : ℕ} (hm : m ≠ 0)
    (x : ↥(powTorsion ((𝒪[K.carrier])ˣ × ZHat) m)) :
    (x : (𝒪[K.carrier])ˣ × ZHat) = ((x : (𝒪[K.carrier])ˣ × ZHat).1, 1) :=
  Prod.ext rfl (torsion_snd_eq_one K hm x)

/-! ## 6. ★作用 —— 大きさだけでは足りない -/

/-- **同型 `e` は `c` 乗作用と可換**。

`𝒪_K^× × Ẑ` の群自己同型 `Φ` が捩れ単数の上で `u ↦ u^c` を誘導するなら、
`μ_m(K)` 側でも `ζ ↦ ζ^c`。★`Φ` が 2 因子を混ぜてもよい(`torsion_eq_mk`)。 -/
theorem torsionUnitsZHatEquiv_pow (K : PAdicLocalField p) {m : ℕ} (hm : m ≠ 0)
    (Φ : ((𝒪[K.carrier])ˣ × ZHat) ≃* ((𝒪[K.carrier])ˣ × ZHat)) (c : ℕ)
    (hΦ : ∀ u : (𝒪[K.carrier])ˣ, u ^ m = 1 → (Φ (u, 1)).1 = u ^ c)
    (x : ↥(powTorsion ((𝒪[K.carrier])ˣ × ZHat) m)) :
    torsionUnitsZHatEquiv K hm (powTorsionCongr Φ m x)
      = (torsionUnitsZHatEquiv K hm x) ^ c := by
  refine Subtype.ext ?_
  have hc : (((torsionUnitsZHatEquiv K hm x) ^ c : ↥(rootsOfUnity m K.carrier)) : (K.carrier)ˣ)
      = (unitsToField K (x : (𝒪[K.carrier])ˣ × ZHat).1) ^ c := by
    rw [SubmonoidClass.coe_pow, torsionUnitsZHatEquiv_coe]
  rw [torsionUnitsZHatEquiv_coe, hc, powTorsionCongr_coe, torsion_eq_mk K hm x,
    hΦ _ (powTorsion_fst x.2), map_pow]

/-- **同型 `e` は Galois 同変**(体の自己同型 `σ` の形)。 -/
theorem torsionUnitsZHatEquiv_galois (K : PAdicLocalField p) {m : ℕ} (hm : m ≠ 0)
    (Φ : ((𝒪[K.carrier])ˣ × ZHat) ≃* ((𝒪[K.carrier])ˣ × ZHat)) (σ : K.carrier ≃+* K.carrier)
    (hΦ : ∀ u : (𝒪[K.carrier])ˣ, u ^ m = 1 →
      ((unitsToField K (Φ (u, 1)).1 : (K.carrier)ˣ) : K.carrier)
        = σ (((u : 𝒪[K.carrier]) : K.carrier)))
    (x : ↥(powTorsion ((𝒪[K.carrier])ˣ × ZHat) m)) :
    (((torsionUnitsZHatEquiv K hm (powTorsionCongr Φ m x) : ↥(rootsOfUnity m K.carrier))
        : (K.carrier)ˣ) : K.carrier)
      = σ ((((torsionUnitsZHatEquiv K hm x : ↥(rootsOfUnity m K.carrier)) : (K.carrier)ˣ))
          : K.carrier) := by
  have h : (Φ (x : (𝒪[K.carrier])ˣ × ZHat)).1
      = (Φ ((x : (𝒪[K.carrier])ˣ × ZHat).1, 1)).1 := by
    rw [← torsion_eq_mk K hm x]
  rw [torsionUnitsZHatEquiv_coe, torsionUnitsZHatEquiv_coe, powTorsionCongr_coe, h]
  exact hΦ _ (powTorsion_fst x.2)

/-- **★`K ⊆ F̄` の中では、`p^n` 乗根への作用は円分指標倍**。

`σ : K ≃ K` が `g ∈ Γ_F` の制限(`ι ∘ σ = g ∘ ι`)なら、`μ_{p^n}(K)` の上で
`σ ζ = ζ^{χ_{F,n}(g)}`。指数の式は `TorsionCyclotomeIsCyclotomic` と同じもの。

証明:`ι` は体からの環準同型なので単射。像で
`cyclotomicCharacter.spec` を使い、引き戻す。 -/
theorem ringEquiv_pow_eq_cyclotomicCharacter (F K : PAdicLocalField p) {n : ℕ}
    (ι : K.carrier →+* F.closure) (g : F.absGal) (σ : K.carrier ≃+* K.carrier)
    (hcompat : ∀ x : K.carrier, ι (σ x) = g (ι x))
    {ζ : K.carrier} (hζ : ζ ^ (p ^ n) = 1) :
    σ ζ = ζ ^ ((PadicInt.toZModPow n
      ((cyclotomicCharacter F.closure p g.toRingEquiv : ℤ_[p]))).val) := by
  refine ι.injective ?_
  rw [hcompat, map_pow]
  exact cyclotomicCharacter.spec (n := n) p g.toRingEquiv (ι ζ)
    (by rw [← map_pow, hζ, map_one])

/-! ## 7. ★到達点(その 2)—— `∃` の形に束ねる -/

theorem pn_ne_zero (p : ℕ) [Fact p.Prime] (n : ℕ) : p ^ n ≠ 0 :=
  pow_ne_zero n (Fact.out : p.Prime).ne_zero

/-- **★★★★★★★★★★★★★★★★★★★★(Λ9)`tors_{p^n}(𝒪_K^× × Ẑ) ≅ μ_{p^n}(K)`、作用つき**。

同型 `e` は `∃` の内側にあり、結論の型は `K` と `n` にしか依存しない
(2026-09-06 の 10 例目の退化——自由なパラメータを結論に出す形——を避けている)。 -/
theorem exists_torsionUnitsZHat_equiv_rootsOfUnity (K : PAdicLocalField p) (n : ℕ) :
    ∃ e : ↥(powTorsion ((𝒪[K.carrier])ˣ × ZHat) (p ^ n)) ≃* ↥(rootsOfUnity (p ^ n) K.carrier),
      ∀ (Φ : ((𝒪[K.carrier])ˣ × ZHat) ≃* ((𝒪[K.carrier])ˣ × ZHat)) (c : ℕ),
        (∀ u : (𝒪[K.carrier])ˣ, u ^ (p ^ n) = 1 → (Φ (u, 1)).1 = u ^ c) →
        ∀ x, e (powTorsionCongr Φ (p ^ n) x) = (e x) ^ c :=
  ⟨torsionUnitsZHatEquiv K (pn_ne_zero p n),
    fun Φ c hΦ x => torsionUnitsZHatEquiv_pow K (pn_ne_zero p n) Φ c hΦ x⟩

/-- **★★★★★★★★★★★★★★★★★★★★★★★★(Λ9、円分指標つき)**。

`K ⊆ F̄` を固定する。`Γ_F` の元 `g` が `σ : K ≃ K` を誘導し、
`𝒪_K^× × Ẑ` の自己同型 `Φ` が捩れの上で `σ` を覆うなら、
`tors_{p^n}(𝒪_K^× × Ẑ) ≅ μ_{p^n}(K)` の下で `Φ` は
**円分指標倍 `ζ ↦ ζ^{χ_{F,n}(g)}`** として働く。

★これが Λ8 の接続点である——Λ8 が局所類体論の相互律から `Φ` と
「`Φ` が `σ` を覆う」を供給すれば、`TorsionCyclotomeIsCyclotomic` が要求する
指数の式がそのまま出る。 -/
theorem exists_torsionUnitsZHat_equiv_cyclotomic (F K : PAdicLocalField p) (n : ℕ)
    (ι : K.carrier →+* F.closure) :
    ∃ e : ↥(powTorsion ((𝒪[K.carrier])ˣ × ZHat) (p ^ n)) ≃* ↥(rootsOfUnity (p ^ n) K.carrier),
      ∀ (g : F.absGal) (σ : K.carrier ≃+* K.carrier)
        (Φ : ((𝒪[K.carrier])ˣ × ZHat) ≃* ((𝒪[K.carrier])ˣ × ZHat)),
        (∀ y : K.carrier, ι (σ y) = g (ι y)) →
        (∀ u : (𝒪[K.carrier])ˣ, u ^ (p ^ n) = 1 →
          ((unitsToField K (Φ (u, 1)).1 : (K.carrier)ˣ) : K.carrier)
            = σ (((u : 𝒪[K.carrier]) : K.carrier))) →
        ∀ x, e (powTorsionCongr Φ (p ^ n) x)
          = (e x) ^ ((PadicInt.toZModPow n
              ((cyclotomicCharacter F.closure p g.toRingEquiv : ℤ_[p]))).val) := by
  refine ⟨torsionUnitsZHatEquiv K (pn_ne_zero p n), fun g σ Φ hcompat hΦ x => ?_⟩
  refine torsionUnitsZHatEquiv_pow K (pn_ne_zero p n) Φ _ (fun u hu => ?_) x
  refine units_eq_pow_of_coe K ?_
  rw [hΦ u hu, map_pow, Units.val_pow_eq_pow_val, unitsToField_coe,
    ringEquiv_pow_eq_cyclotomicCharacter F K ι g σ hcompat
      (coe_pow_eq_one_of_units_pow_eq_one K hu)]

end ABC3.Found.PGC
