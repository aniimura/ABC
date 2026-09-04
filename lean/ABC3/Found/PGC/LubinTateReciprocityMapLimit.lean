import ABC3.Found.PGC.LubinTateReciprocityLimitCompat
import ABC3.Found.PGC.LubinTateReciprocityHomomorphism
import ABC3.Found.PGC.UnitsInverseLimit

/-!
# `reciprocityMapLimit : Gal(K.closure/K.carrier) →* CompatibleUnits K hπmax`(`sorry` 無し)

節目(5)(射影極限 `Gal(L_π/K)≅𝒪_K^×`)の最終組み立て、部品(ii)。
`Found/PGC/LubinTateReciprocityLimitCompat.lean::reciprocityMapLimit
Compat`((𝒪_K)^×⧸principalUnits レベルでのn跨ぎ両立性)・
`Found/PGC/LubinTateGeneratorSequence.lean::psiGenSeq`(無限compatible
列)・`Found/PGC/UnitsInverseLimit.lean::principalUnitsQuotientEquiv_
apply_mk`/`principalUnitsQuotientEquiv_map_eq`(2つの自然性)・
`Found/PGC/LubinTateReciprocityHomomorphism.lean::reciprocityMap_mul`
(準同型性)——すべてsorry無しで既に確立済みの部品——を実際に
組み合わせ、`Gal(K.closure/K.carrier)→CompatibleUnits K hπmax`という
**群準同型**を構成する。

## 構成の骨格

1. `reciprocityMapLimitFamily σ n`: `n=0`では(`(𝒪_K/π^0)^×`が自明群
   なので)`1`、`n=m+1`では`psiGenSeq m`(level`m+1`の生成元)の上で
   評価した`reciprocityMap`を`principalUnitsQuotientEquiv`で`(𝒪_K/
   π^{m+1})^×`へ落としたもの。
2. `reciprocityMapLimitFamily_step`: 隣接するn(`n`と`n+1`)の間の
   両立性——`n=0`は自明群どうしなので`Subsingleton.elim`、`n=m+1`は
   `principalUnitsQuotientEquiv_map_eq`(遷移写像との可換性)と
   `reciprocityMapLimitCompat`(n跨ぎ両立性の核心)を組み合わせる。
   ここで唯一の技術的な障害は、`reciprocityMapLimitCompat`が
   `psiGenStepResult m (psiGenSeq m)`という形で生成元を参照するのに
   対し、`reciprocityMapLimitFamily`は`psiGenSeq (m+1)`という形で
   参照すること——両者は`rfl`で一致する(`psiGenSeq`の定義から)が、
   `rw`は**構文的一致**しか見ないため素通りしてしまう。★対処:
   `rw`ではなく`congrArg`(defeqチェックで通る)を使う。
3. `compatible_of_succ`(完全に一般的な、ℕ添字付き両立系に関する
   補題——`α`・遷移写像`tr`が満たすべき合成規則(`htrans`)・恒等規則
   (`hrefl`)・隣接両立性(`hstep`)から、**任意の`m≤n`**での両立性を
   `Nat.le_induction`で導く): これを`unitReductionTransition`
   (合成規則`unitReductionTransition_trans`・恒等規則
   `unitReductionTransition_refl`、ともに`Ideal.Quotient.factor_comp_
   apply`/`factor_mk`から直ちに従う)に適用し、`reciprocityMapLimit
   Family_step`(隣接両立性)から`reciprocityMapLimitFamily_compatible`
   (任意の`m≤n`での両立性)を得る——これで`reciprocityMapLimitFamily σ`
   が`CompatibleUnits K hπmax`の元になる。
4. `reciprocityMapLimitFamily_one`/`_mul`: `principalUnitsQuotientEquiv`
   が`MulEquiv`(`map_one`/`map_mul`を自動的に持つ)であることと
   `reciprocityMap_one`/`reciprocityMap_mul`(既出、本セッション以前の
   集大成)を組み合わせるだけ。
5. `reciprocityMapLimit`: 上記をすべて束ね、`toFun`・`map_one'`・
   `map_mul'`を与える`MonoidHom`として組み立てる。

## 逸脱(記録)

`reciprocityMapLimitFamily`の`n=0`成分は数学的には「射影系の最初の
座標(自明群への写像)」であり、原典では暗黙のうちに`n≥1`だけを
扱っている。`CompatibleUnits`の台が`∀n:ℕ,...`(`n=0`も含む)である
以上、形式化ではこの成分を明示的に`1`と定める必要があった——
`(𝒪_K/π^0)^×`が自明群であることから両立性は自動的に満たされ、
後続の証明(全射性・単射性)に影響しないので、原典の意図を優先する
軽微な逸脱として扱う。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped NormedField Valued Classical

/-- **完全に一般的な事実**: ℕ添字付きの族`α`上の遷移写像`tr`が合成
規則(`htrans`)・恒等規則(`hrefl`)を満たし、かつ**隣接する**添字の
間で両立系`v`と両立する(`hstep`)なら、**任意の`m≤n`**についても
両立する。`Nat.le_induction`(`m`を固定して`n`を`m`から上に動かす
帰納法)で直ちに従う——`unitReductionTransition`に限らず再利用できる
形にしてある。 -/
theorem compatible_of_succ {α : ℕ → Type*} (tr : ∀ {m n : ℕ}, m ≤ n → α n → α m)
    (htrans : ∀ {m n k : ℕ} (h1 : m ≤ n) (h2 : n ≤ k) (x : α k), tr h1 (tr h2 x) = tr (h1.trans h2) x)
    (hrefl : ∀ {n : ℕ} (x : α n), tr (le_refl n) x = x)
    (v : ∀ n, α n) (hstep : ∀ n, tr (Nat.le_succ n) (v (n + 1)) = v n) :
    ∀ {m n : ℕ} (h : m ≤ n), tr h (v n) = v m := by
  intro m n h
  induction n, h using Nat.le_induction with
  | base => exact hrefl (v m)
  | succ n hmn ih =>
    rw [← htrans hmn (Nat.le_succ n) (v (n + 1)), hstep n, ih]

/-- `unitReductionTransition`の恒等規則: `h=le_refl n`のときは恒等
写像——`Ideal.Quotient.factor_mk`から`mk`の像の上で確認すれば十分。 -/
theorem unitReductionTransition_refl {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (n : ℕ) (x : (𝒪[K.carrier] ⧸ Ideal.span ({π ^ n} : Set (𝒪[K.carrier])))ˣ) :
    unitReductionTransition K hπmax (le_refl n) x = x := by
  apply Units.ext
  unfold unitReductionTransition
  obtain ⟨a, ha⟩ := Ideal.Quotient.mk_surjective (x : 𝒪[K.carrier] ⧸ Ideal.span ({π ^ n} : Set (𝒪[K.carrier])))
  show Ideal.Quotient.factor (le_refl _) (x : 𝒪[K.carrier] ⧸ Ideal.span ({π ^ n} : Set (𝒪[K.carrier]))) = x
  rw [← ha, Ideal.Quotient.factor_mk]

/-- `unitReductionTransition`の合成規則: `m≤n≤k`の2段の遷移は
`m≤k`への直接の遷移に一致する——`Ideal.Quotient.factor_mk`を`mk`の
像の上で3回使うだけ(`factor`同士の合成は`Ideal.Quotient.factor_comp`
から従うが、ここでは`mk`経由の等式で十分)。 -/
theorem unitReductionTransition_trans {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    {m n k : ℕ} (h1 : m ≤ n) (h2 : n ≤ k)
    (x : (𝒪[K.carrier] ⧸ Ideal.span ({π ^ k} : Set (𝒪[K.carrier])))ˣ) :
    unitReductionTransition K hπmax h1 (unitReductionTransition K hπmax h2 x) =
      unitReductionTransition K hπmax (h1.trans h2) x := by
  apply Units.ext
  unfold unitReductionTransition
  obtain ⟨a, ha⟩ := Ideal.Quotient.mk_surjective (x : 𝒪[K.carrier] ⧸ Ideal.span ({π ^ k} : Set (𝒪[K.carrier])))
  show Ideal.Quotient.factor
      (S := Ideal.span ({π ^ n} : Set (𝒪[K.carrier]))) (T := Ideal.span ({π ^ m} : Set (𝒪[K.carrier])))
      (by rw [Ideal.span_singleton_le_span_singleton]; exact pow_dvd_pow π h1)
      (Ideal.Quotient.factor
        (S := Ideal.span ({π ^ k} : Set (𝒪[K.carrier]))) (T := Ideal.span ({π ^ n} : Set (𝒪[K.carrier])))
        (by rw [Ideal.span_singleton_le_span_singleton]; exact pow_dvd_pow π h2)
        (x : 𝒪[K.carrier] ⧸ Ideal.span ({π ^ k} : Set (𝒪[K.carrier])))) =
    Ideal.Quotient.factor
      (S := Ideal.span ({π ^ k} : Set (𝒪[K.carrier]))) (T := Ideal.span ({π ^ m} : Set (𝒪[K.carrier])))
      (by rw [Ideal.span_singleton_le_span_singleton]; exact pow_dvd_pow π (h1.trans h2))
      (x : 𝒪[K.carrier] ⧸ Ideal.span ({π ^ k} : Set (𝒪[K.carrier])))
  rw [← ha, Ideal.Quotient.factor_mk, Ideal.Quotient.factor_mk, Ideal.Quotient.factor_mk]

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**`reciprocityMapLimit`本体の族**: 各`σ`・各レベル`n`について
`(𝒪_K/π^n)^×`の元を与える。`n=0`は自明群への`1`(逸脱として上の
module docstringに記録)、`n=m+1`は`psiGenSeq m`(level`m+1`の
compatible生成元)の上で評価した`reciprocityMap`を
`principalUnitsQuotientEquiv`で落としたもの。 -/
noncomputable def reciprocityMapLimitFamily
    {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (σ : K.closure ≃ₐ[K.carrier] K.closure) :
    (n : ℕ) → (𝒪[K.carrier] ⧸ Ideal.span ({π ^ n} : Set (𝒪[K.carrier])))ˣ
  | 0 => 1
  | (n + 1) =>
    haveI := (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf n).hfd
    principalUnitsQuotientEquiv K hπmax (n + 1) (by omega)
      (reciprocityMap K hq hπmax hπne0 f hf0 hf1 hf (n + 1) (by omega)
        (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf n).pt
        (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf n).hψ
        (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf n).hn
        (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf n).hmem σ)

/-- **`reciprocityMapLimitFamily`の隣接両立性**: `n`と`n+1`の間で
`unitReductionTransition`と可換。`n=0`は`(𝒪_K/π^0)^×`が自明群である
ことから、`n=m+1`は`principalUnitsQuotientEquiv_map_eq`+
`reciprocityMapLimitCompat`の組み合わせから従う——`congrArg`で
`psiGenStepResult`基準と`psiGenSeq`基準の(`rfl`で一致する)表現の
差を吸収する。 -/
theorem reciprocityMapLimitFamily_step
    {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (σ : K.closure ≃ₐ[K.carrier] K.closure) (n : ℕ) :
    unitReductionTransition K hπmax (Nat.le_succ n)
        (reciprocityMapLimitFamily K hq hπmax hπne0 f hf0 hf1 hf σ (n + 1)) =
      reciprocityMapLimitFamily K hq hπmax hπne0 f hf0 hf1 hf σ n := by
  match n with
  | 0 =>
    have h : Subsingleton (𝒪[K.carrier] ⧸ Ideal.span ({π ^ 0} : Set (𝒪[K.carrier])))ˣ := by
      have h0 : Ideal.span ({π ^ 0} : Set (𝒪[K.carrier])) = ⊤ := by
        rw [pow_zero]; exact Ideal.span_singleton_one
      rw [h0]
      constructor
      intro a b
      apply Units.ext
      exact Subsingleton.elim _ _
    exact h.elim _ _
  | (m + 1) =>
    haveI := (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf (m + 1)).hfd
    haveI := (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hfd
    haveI := (psiGenStepResult K hq hπmax hπne0 f hf0 hf1 hf m (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m)).hfd
    show unitReductionTransition K hπmax (Nat.le_succ (m + 1))
        (principalUnitsQuotientEquiv K hπmax (m + 1 + 1) (by omega)
            (reciprocityMap K hq hπmax hπne0 f hf0 hf1 hf (m + 1 + 1) (by omega)
              (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf (m + 1)).pt
              (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf (m + 1)).hψ
              (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf (m + 1)).hn
              (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf (m + 1)).hmem σ)) =
      principalUnitsQuotientEquiv K hπmax (m + 1) (by omega)
          (reciprocityMap K hq hπmax hπne0 f hf0 hf1 hf (m + 1) (by omega)
            (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt
            (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hψ
            (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hn
            (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hmem σ)
    rw [← principalUnitsQuotientEquiv_map_eq K hπmax (m + 1) (by omega)]
    exact congrArg (principalUnitsQuotientEquiv K hπmax (m + 1) (by omega))
      (reciprocityMapLimitCompat K hq hπmax hπne0 f hf0 hf1 hf m σ)

/-- **`reciprocityMapLimitFamily`は`CompatibleUnits`の元を与える**:
任意の`m≤n`について両立する。`compatible_of_succ`を
`unitReductionTransition`(合成規則・恒等規則を上の2補題で確認済み)
に適用し、隣接両立性(`reciprocityMapLimitFamily_step`)から得る。 -/
theorem reciprocityMapLimitFamily_compatible
    {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (σ : K.closure ≃ₐ[K.carrier] K.closure) {m n : ℕ} (h : m ≤ n) :
    unitReductionTransition K hπmax h (reciprocityMapLimitFamily K hq hπmax hπne0 f hf0 hf1 hf σ n) =
      reciprocityMapLimitFamily K hq hπmax hπne0 f hf0 hf1 hf σ m :=
  compatible_of_succ (fun {_ _} h x => unitReductionTransition K hπmax h x)
    (fun h1 h2 x => unitReductionTransition_trans K hπmax h1 h2 x)
    (fun x => unitReductionTransition_refl K hπmax _ x)
    (reciprocityMapLimitFamily K hq hπmax hπne0 f hf0 hf1 hf σ)
    (reciprocityMapLimitFamily_step K hq hπmax hπne0 f hf0 hf1 hf σ) h

/-- **`reciprocityMapLimitFamily`は`σ=1`で恒等的に`1`**——
`reciprocityMap_one`(既出)と`principalUnitsQuotientEquiv`が
`MulEquiv`(`map_one`を自動的に持つ)であることから。 -/
theorem reciprocityMapLimitFamily_one
    {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) :
    reciprocityMapLimitFamily K hq hπmax hπne0 f hf0 hf1 hf 1 n = 1 := by
  match n with
  | 0 => rfl
  | (m + 1) =>
    haveI := (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hfd
    show principalUnitsQuotientEquiv K hπmax (m + 1) (by omega)
        (reciprocityMap K hq hπmax hπne0 f hf0 hf1 hf (m + 1) (by omega)
          (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt
          (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hψ
          (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hn
          (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hmem 1) = 1
    rw [show (1 : K.closure ≃ₐ[K.carrier] K.closure) = AlgEquiv.refl from rfl,
      reciprocityMap_one]
    exact map_one _

/-- **`reciprocityMapLimitFamily`は準同型**: `σ*τ`での値は`σ`・`τ`
それぞれの値の積。`reciprocityMap_mul`(既出、本セッション以前の
集大成)と`principalUnitsQuotientEquiv`が`MulEquiv`であること
(`map_mul`)から。 -/
theorem reciprocityMapLimitFamily_mul
    {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (σ τ : K.closure ≃ₐ[K.carrier] K.closure) (n : ℕ) :
    reciprocityMapLimitFamily K hq hπmax hπne0 f hf0 hf1 hf (σ * τ) n =
      reciprocityMapLimitFamily K hq hπmax hπne0 f hf0 hf1 hf σ n *
        reciprocityMapLimitFamily K hq hπmax hπne0 f hf0 hf1 hf τ n := by
  match n with
  | 0 => exact (mul_one 1).symm
  | (m + 1) =>
    haveI := (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hfd
    show principalUnitsQuotientEquiv K hπmax (m + 1) (by omega)
        (reciprocityMap K hq hπmax hπne0 f hf0 hf1 hf (m + 1) (by omega)
          (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt
          (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hψ
          (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hn
          (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hmem (σ * τ)) =
      principalUnitsQuotientEquiv K hπmax (m + 1) (by omega)
          (reciprocityMap K hq hπmax hπne0 f hf0 hf1 hf (m + 1) (by omega)
            (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt
            (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hψ
            (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hn
            (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hmem σ) *
        principalUnitsQuotientEquiv K hπmax (m + 1) (by omega)
          (reciprocityMap K hq hπmax hπne0 f hf0 hf1 hf (m + 1) (by omega)
            (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt
            (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hψ
            (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hn
            (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hmem τ)
    rw [reciprocityMap_mul, map_mul]

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**`reciprocityMapLimit`**: `Gal(K.closure/K.carrier)→CompatibleUnits
K hπmax`という**群準同型**。節目(5)(射影極限`Gal(L_π/K)≅𝒪_K^×`)の
最終組み立て、部品(ii)の完成——`reciprocityMapLimitFamily`が各`σ`に
ついて`CompatibleUnits`の元を与えること(`_compatible`)と、準同型性
(`_one`・`_mul`)を束ねるだけ。残るのは(iii)全射性・単射性のみ。 -/
noncomputable def reciprocityMapLimit
    {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff)) :
    (K.closure ≃ₐ[K.carrier] K.closure) →* CompatibleUnits K hπmax where
  toFun σ := ⟨reciprocityMapLimitFamily K hq hπmax hπne0 f hf0 hf1 hf σ,
    fun {m n} h => reciprocityMapLimitFamily_compatible K hq hπmax hπne0 f hf0 hf1 hf σ h⟩
  map_one' := by
    apply Subtype.ext
    funext n
    exact reciprocityMapLimitFamily_one K hq hπmax hπne0 f hf0 hf1 hf n
  map_mul' := by
    intro σ τ
    apply Subtype.ext
    funext n
    exact reciprocityMapLimitFamily_mul K hq hπmax hπne0 f hf0 hf1 hf σ τ n

end ABC3.Found.PGC
