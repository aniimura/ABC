import ABC3.Found.PGC.KummerDuality
import ABC3.Found.PGC.UnramifiedClosureRoots
import ABC3.Found.PGC.AdjoinFieldClosure
import ABC3.Found.PGC.ResidueCardLowerBound
import ABC3.Found.PGC.InertiaIdentification

/-!
# 惰性群からの連続準同型 —— (C-q) の上界

[pGC] Proposition 1.2 への経路 C(`ResearchPaper/pgc-goal.md` のノード F)。
`N_n(Γ_K) = # Hom_cont(Γ_K, ℤ/n)` の**上界**を出すために、惰性群

  `I_K = Gal(K̄/K^ur) = (unramifiedClosure K).fixingSubgroup`

からの連続準同型を Kummer 理論で調べる。本ファイルの到達点は

  `pow_residueCard_sub_one_restrict` ——
  `f : Γ_K →* ℤ/n` が連続(核が開)で `p ∤ n` なら、
  `σ ∈ I_K` に対して `f σ ^ (q − 1) = 1`。

すなわち **`f` の惰性部分は `gcd(n, q−1)` 次の部分群にしか値を持てない**。

もう一つの到達点は

  `contHom_fixingSubgroup_unramifiedClosure`(F1) ——
  `p ∤ n` なら `Hom_cont(I_K, ℤ/n)` は**位数 `n` の巡回群**。

両者を合わせて `card_range_restrictInertia_le`
(惰性群への制限射の像は `gcd(n, q−1)` 個以下)が出る。

★**本ファイルに無いもの(F3)**:`Hom_cont(Gal(K^ur/K), ℤ/n)` の個数が `n` であること。
(C-q) の上界 `N_n(Γ_K) ≤ n·gcd(n, q−1)` はそれだけを仮定に残した形
`contHomCard_absGal_le_of_card_ker` で用意してある。

## ★惰性群を `inertia` で書かない(設計上の要点)

`Found/PGC/InertiaIdentification.lean` の `inertia` は `Interface` の
`SubgroupCorrespondence` / `ResidueCardinality` を経由する。経路 C の出口は
Corollary 1.3(惰性群が `Γ_K` から決まる)なので、そこを通ると**循環する**。
本ファイルは惰性群を常に `(unramifiedClosure K).fixingSubgroup` と直接書く。
`inertia` も `inertia_eq_fixingSubgroup_unramifiedClosure` も参照しない。

## 証明の骨格(`pow_residueCard_sub_one_restrict`)

1. **Kummer(全射性)**——`μ_n ⊆ K^ur`(在庫、第 1005/1006)なので、
   `F := K^ur` 上の Kummer 双対(`Found/PGC/KummerDuality.lean`、第 1020)が使える。
   `f` を `I_K` に制限し `Gal(K̄/K^ur)` の指標と見ると、ある `β ∈ K̄` が取れて
   `β^n ∈ K^ur` かつ **`f` は `σ ↦ σβ/β` を通してしか `I_K` を見ない**
   (`exists_kummer_root_inertia`)。
2. **`τβ/β` は `K^ur` の中にある**——`b := β^n` とすると `‖τ b / b‖ = 1`
   (ノルムは Galois 不変)なので、在庫
   `exists_pow_eq_mem_unramifiedClosure`(第 1018、`p ∤ n` なら `K^ur` の
   ノルム 1 の元は `K^ur` の中で `n` 乗数)から `τ b / b = w^n`(`w ∈ K^ur`)。
   `v := τβ/β` は `(v/w)^n = 1` をみたすので `v/w ∈ μ_n ⊆ K^ur`、よって
   **`v ∈ K^ur`**。★ここが要で、`K^ur` の付値群を調べずに済ませている。
3. **共役**——`σ' := τ⁻¹στ` とすると `σ'β/β = τ⁻¹(σβ/β)`(2 の `σ v = v` から)。
4. **Frobenius は `μ_n` に `q` 乗で作用**(在庫、第 1010)ので、`κ := σβ/β` に対し
   `κ = (τ⁻¹κ)^q`。したがって `σ'^q` と `σ` は 1 の同じ Kummer 値を与え、
   1 から `f (σ'^q) = f σ`。左辺は `f σ ^ q`(`f` は可換群に値を取るので
   共役で不変)なので `f σ ^ q = f σ`、すなわち `f σ ^ (q−1) = 1`。

★原典(および Serre の局所類体論)は相互律を経由するが、本経路は経由しない。
これは `pgc-goal.md` に記録済みの逸脱である。

## ★逸脱の記録

* 惰性群を `inertia` ではなく `(unramifiedClosure K).fixingSubgroup` と書く(上記)。
* 「連続」を `IsOpen (ker f)` で表す(`Found/PGC/ContinuousHomCount.lean` の規約。
  `ZMod n` に `TopologicalSpace` インスタンスが無いための措置)。
-/

namespace ABC3.Found.PGC

/-! ## Kummer 全射性を「根の形」で書き直し、任意の代数閉包へ移す

`Found/PGC/KummerDuality.lean` の `range_kummerMap` は係数群 `μ_n(F̄)` と
**`AlgebraicClosure F` という特定の型**で書かれている。消費側では
`Ω = K.closure`(`F = K^ur` 上の代数閉包でもある)に当てたいので、
ここで 2 段の書き換えをする:

1. 係数群を `Multiplicative (ZMod n)` に替え、指標を「`α` の動き」だけで書く。
2. `AlgEquiv.autCongr`(第 993 で連続性が済んでいる)で任意の代数閉包へ移す。
-/

section KummerRoot

open ABC3.Found.PGC.KummerDual

variable {F : Type} [Field F] [CharZero F] {n : ℕ} {ζ : F}

/-- **Kummer 全射性(根の形、`AlgebraicClosure F` 版)**。

`μ_n ⊆ F` のとき、`Γ_F` の連続指標 `g : Γ_F →* ℤ/n` はある `α`(`α^n ∈ F`)を
使って `σ ↦ σα/α` を通してしか `Γ_F` を見ない。

★係数群を `μ_n` でなく `ℤ/n` のまま扱い、指標そのものではなく
「`σα/α` が等しければ値が等しい」という**弱い形**にしてある。これで
`μ_n(Ω) ≃* ℤ/n` の同一視を消費側に持ち込まずに済む。

退化の自己検査:`hn : n ≠ 0` を落とすと `ZMod 0 = ℤ` で `μ_0` が定まらず
`hζ` が成り立たない(`IsPrimitiveRoot ζ 0` は `ζ` が単元であることしか言わない)ので
主張は偽になる。`hζ`(原始根が `F` にある)を落とすと Kummer 理論そのものが
壊れる——`F = ℚ`, `n = 3` では `Γ_ℚ → ℤ/3` の連続指標(例えば `ℚ(ζ_9)^+` の
巡回 3 次拡大が与えるもの)が `σα/α` の形に書けない。`CharZero F` は
`KummerDuality` 側の要求をそのまま引き継いだもの。 -/
theorem exists_root_of_contHom_algClosure (hζ : IsPrimitiveRoot ζ n) (hn : n ≠ 0)
    (g : (AlgebraicClosure F ≃ₐ[F] AlgebraicClosure F) →* Multiplicative (ZMod n))
    (hg : IsOpen ((MonoidHom.ker g : Subgroup (AlgebraicClosure F ≃ₐ[F] AlgebraicClosure F)) :
      Set (AlgebraicClosure F ≃ₐ[F] AlgebraicClosure F))) :
    ∃ α : AlgebraicClosure F, α ≠ 0 ∧
      (∃ c : F, α ^ n = algebraMap F (AlgebraicClosure F) c) ∧
      ∀ σ τ : AlgebraicClosure F ≃ₐ[F] AlgebraicClosure F, σ α / α = τ α / α → g σ = g τ := by
  set e := zmodMulEquivRootsOfUnity (F := F) hζ hn with he
  set h : (AlgebraicClosure F ≃ₐ[F] AlgebraicClosure F) →*
      ↥(rootsOfUnity n (AlgebraicClosure F)) := e.toMonoidHom.comp g with hh
  have hmem : h ∈ contHom (AlgebraicClosure F ≃ₐ[F] AlgebraicClosure F)
      ↥(rootsOfUnity n (AlgebraicClosure F)) := by
    rw [mem_contHom, hh, MonoidHom.ker_comp_of_injective _ _ e.injective]
    exact hg
  rw [← range_kummerMap hζ hn] at hmem
  obtain ⟨a, ha⟩ := hmem
  refine ⟨nthRoot hn a, ne_zero_of_pow_eq hn (nthRoot_pow hn a), ⟨(a : F), nthRoot_pow hn a⟩,
    fun σ τ hst => ?_⟩
  have key : h σ = h τ := by
    rw [← ha]
    apply Subtype.ext
    apply Units.ext
    show σ (nthRoot hn a) / nthRoot hn a = τ (nthRoot hn a) / nthRoot hn a
    exact hst
  exact e.injective (by simpa [hh] using key)

/-- **Kummer 全射性(根の形、任意の代数閉包 `Ω` 版)**。

`AlgEquiv.autCongr`(`Found/PGC/AdjoinFieldClosure.lean`、第 993 で連続性済み)で
`AlgebraicClosure F` から `Ω` へ移しただけ。消費側(`Ω = K.closure`)では
`AlgebraicClosure ↥(K^ur)` という**別の型**を一切扱わずに済む。

退化の自己検査:`[IsAlgClosure F Ω]` を落とすと `Ω` を `F` の適当な有限次拡大に
取れてしまい、`Gal(Ω/F)` が小さすぎて偽になる。他の仮定は上の
`exists_root_of_contHom_algClosure` と同じ。 -/
theorem exists_root_of_contHom {Ω : Type} [Field Ω] [Algebra F Ω] [IsAlgClosure F Ω]
    (hζ : IsPrimitiveRoot ζ n) (hn : n ≠ 0)
    (g : (Ω ≃ₐ[F] Ω) →* Multiplicative (ZMod n))
    (hg : IsOpen ((MonoidHom.ker g : Subgroup (Ω ≃ₐ[F] Ω)) : Set (Ω ≃ₐ[F] Ω))) :
    ∃ β : Ω, β ≠ 0 ∧ (∃ c : F, β ^ n = algebraMap F Ω c) ∧
      ∀ σ τ : Ω ≃ₐ[F] Ω, σ β / β = τ β / β → g σ = g τ := by
  set e : AlgebraicClosure F ≃ₐ[F] Ω := IsAlgClosure.equiv F _ _ with he
  set Φ : (AlgebraicClosure F ≃ₐ[F] AlgebraicClosure F) ≃* (Ω ≃ₐ[F] Ω) :=
    AlgEquiv.autCongr e with hΦ
  have hΦapp : ∀ (ψ : AlgebraicClosure F ≃ₐ[F] AlgebraicClosure F) (y : AlgebraicClosure F),
      (Φ ψ) (e y) = e (ψ y) := by
    intro ψ y
    show e (ψ (e.symm (e y))) = e (ψ y)
    simp
  set g' : (AlgebraicClosure F ≃ₐ[F] AlgebraicClosure F) →* Multiplicative (ZMod n) :=
    g.comp Φ.toMonoidHom with hg'
  have hg'open : IsOpen ((MonoidHom.ker g' :
      Subgroup (AlgebraicClosure F ≃ₐ[F] AlgebraicClosure F)) :
      Set (AlgebraicClosure F ≃ₐ[F] AlgebraicClosure F)) :=
    comp_mem_contHom (A := Multiplicative (ZMod n)) (Φ.toMonoidHom)
      (continuous_autCongr e) (f := g) hg
  obtain ⟨α, hα0, ⟨c, hc⟩, hkey⟩ := exists_root_of_contHom_algClosure hζ hn g' hg'open
  refine ⟨e α, by simpa using hα0, ⟨c, by rw [← map_pow, hc, AlgEquiv.commutes]⟩, ?_⟩
  intro σ τ hst
  have h1 : (Φ.symm σ) α / α = (Φ.symm τ) α / α := by
    apply e.injective
    rw [map_div₀, map_div₀, ← hΦapp, ← hΦapp, Φ.apply_symm_apply, Φ.apply_symm_apply]
    exact hst
  have h2 := hkey _ _ h1
  simpa [hg'] using h2

end KummerRoot

/-! ## `K^ur` 上の Kummer 理論 -/

open ABC3.Skeleton.PGC
open scoped NormedField Valued

variable {p : ℕ} [Fact p.Prime]

/-- `K.closure` は `K^ur` 上の**代数閉包**でもある
(`AdjoinFieldClosure.lean` の `isAlgClosureAdjoinField` と同じ形)。 -/
instance isAlgClosureUnramifiedClosure (K : PAdicLocalField p) :
    IsAlgClosure ↥(unramifiedClosure K) K.closure := ⟨inferInstance, inferInstance⟩

/-- **★★★★★★★★連続指標は惰性群を「`σβ/β`」を通してしか見ない**。

`p ∤ n` なら `μ_n ⊆ K^ur`(在庫)なので `F = K^ur` 上で Kummer 双対が使える。
`f : Γ_K →* ℤ/n` を `I_K` に制限して `Gal(K̄/K^ur)` の指標と見れば、
`exists_root_of_contHom` がその指標を `σ ↦ σβ/β` の形に書き下す。

退化の自己検査:`hn : ¬ p ∣ n` を落とすと `μ_n ⊄ K^ur`(`μ_p` を含む拡大は
分岐する)ので Kummer 双対の前提が崩れ、結論は偽になる——`p ∣ n` のときの
`Hom_cont(I_K, ℤ/n)` は野性的惰性群の寄与で巨大であり、`K^ur` の一つの元
`β^n` では書き尽くせない。`hf`(核が開)を落とすと `range_kummerMap` の
上界側が使えず、証明は届かない。 -/
theorem exists_kummer_root_inertia (K : PAdicLocalField p) {n : ℕ} (hn : ¬ p ∣ n)
    (f : K.absGal →* Multiplicative (ZMod n))
    (hf : IsOpen ((MonoidHom.ker f : Subgroup K.absGal) : Set K.absGal)) :
    ∃ β : K.closure, β ≠ 0 ∧ β ^ n ∈ unramifiedClosure K ∧
      ∀ σ τ : K.absGal, σ ∈ (unramifiedClosure K).fixingSubgroup →
        τ ∈ (unramifiedClosure K).fixingSubgroup → σ β / β = τ β / β → f σ = f τ := by
  have hn0 : n ≠ 0 := by rintro rfl; exact hn (dvd_zero p)
  set E := unramifiedClosure K with hE
  obtain ⟨ζ0, hmem, hζ0⟩ := exists_isPrimitiveRoot_mem_unramifiedClosure K n
    (Nat.pos_of_ne_zero hn0) hn
  have hζ : IsPrimitiveRoot (⟨ζ0, hmem⟩ : ↥E) n :=
    IsPrimitiveRoot.of_map_of_injective (f := algebraMap ↥E K.closure) hζ0
      (algebraMap ↥E K.closure).injective
  have hres : f.comp (E.fixingSubgroup).subtype
      ∈ contHom ↥(E.fixingSubgroup) (Multiplicative (ZMod n)) := by
    rw [mem_contHom]
    have h : ((MonoidHom.ker (f.comp (E.fixingSubgroup).subtype) :
        Subgroup ↥(E.fixingSubgroup)) : Set ↥(E.fixingSubgroup))
        = Subtype.val ⁻¹' ((MonoidHom.ker f : Subgroup K.absGal) : Set K.absGal) := by
      ext σ; simp [MonoidHom.mem_ker]
    rw [h]
    exact hf.preimage continuous_subtype_val
  have hgmem : ((f.comp (E.fixingSubgroup).subtype).comp
      E.fixingSubgroupEquiv.symm.toMonoidHom)
      ∈ contHom (K.closure ≃ₐ[↥E] K.closure) (Multiplicative (ZMod n)) :=
    comp_mem_contHom (A := Multiplicative (ZMod n)) _
      (continuous_fixingSubgroupEquivInf_symm E) hres
  rw [mem_contHom] at hgmem
  obtain ⟨β, hβ0, ⟨c, hc⟩, hkey⟩ := exists_root_of_contHom hζ hn0 _ hgmem
  refine ⟨β, hβ0, ?_, ?_⟩
  · rw [hc]; exact c.2
  · intro σ τ hσ hτ hst
    have h1 := hkey (E.fixingSubgroupEquiv ⟨σ, hσ⟩) (E.fixingSubgroupEquiv ⟨τ, hτ⟩) hst
    simpa using h1

/-- **★★★★★★`μ_n` に `q` 乗として働く `Γ_K` の元(Frobenius の持ち上げ)**。

原始 `n` 乗根 `ζ` は `K^ur` にある(在庫)ので、ある不分岐 `K(x)/K` に含まれる。
その次数 `d` の不分岐拡大は一意(在庫 `adjoin_eq_of_isUnramified`)なので、
`exists_frobenius_pow_of_pow_eq_one`(在庫、第 1010)が与える次数 `d` の
Frobenius `σ` は `ζ` に `q` 乗で働く。`K(x)/K` は normal なので
`AlgEquiv.restrictNormalHom` の全射性で `Γ_K` へ持ち上がる。

退化の自己検査:`hn : ¬ p ∣ n` を落とすと `μ_n ⊄ K^ur` になり、`ζ` を不分岐拡大に
取れないので構成が壊れる。実際 `n = p` では `μ_p ⊆ K(μ_p)` は完全分岐であり、
`Γ_K` の元が `μ_p` に `q` 乗で働くとは限らない(円分指標が非自明)。 -/
theorem exists_absGal_pow_residueCard (K : PAdicLocalField p) {n : ℕ} (hn : ¬ p ∣ n) :
    ∃ τ : K.absGal, ∀ z : K.closure, z ^ n = 1 → τ z = z ^ (Nat.card 𝓀[K.carrier]) := by
  have hn0 : n ≠ 0 := by rintro rfl; exact hn (dvd_zero p)
  haveI : NeZero n := ⟨hn0⟩
  obtain ⟨ζ0, hmem, hζ0⟩ := exists_isPrimitiveRoot_mem_unramifiedClosure K n
    (Nat.pos_of_ne_zero hn0) hn
  obtain ⟨x, hux, hζmem⟩ := (mem_unramifiedClosure_iff K ζ0).mp hmem
  set d := Module.finrank K.carrier
    (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) with hd
  have hd0 : d ≠ 0 := Module.finrank_pos.ne'
  obtain ⟨y, hrank, huy, σ, -, -, hpow⟩ := exists_frobenius_pow_of_pow_eq_one K d hd0
  have heq : IntermediateField.adjoin K.carrier ({x} : Set K.closure)
      = IntermediateField.adjoin K.carrier ({y} : Set K.closure) :=
    adjoin_eq_of_isUnramified K x y hux huy (by rw [hrank])
  have hζy : ζ0 ∈ IntermediateField.adjoin K.carrier ({y} : Set K.closure) := heq ▸ hζmem
  haveI : Normal K.carrier (IntermediateField.adjoin K.carrier ({y} : Set K.closure)) :=
    normal_of_isUnramifiedAdjoin K y huy
  haveI : Normal K.carrier K.closure := IsAlgClosure.normal K.carrier K.closure
  obtain ⟨τ, hτ⟩ := AlgEquiv.restrictNormalHom_surjective
    (F := K.carrier) (K₁ := IntermediateField.adjoin K.carrier ({y} : Set K.closure))
    (E := K.closure) σ
  have hτζ : τ ζ0 = ζ0 ^ (Nat.card 𝓀[K.carrier]) := by
    have h1 : σ ⟨ζ0, hζy⟩ = (⟨ζ0, hζy⟩ :
        IntermediateField.adjoin K.carrier ({y} : Set K.closure)) ^ (Nat.card 𝓀[K.carrier]) := by
      refine hpow hn _ ?_
      exact Subtype.ext (by simpa using hζ0.pow_eq_one)
    have h2 := AlgEquiv.restrictNormal_commutes τ
      (IntermediateField.adjoin K.carrier ({y} : Set K.closure)) ⟨ζ0, hζy⟩
    have h3 : AlgEquiv.restrictNormal τ
        (IntermediateField.adjoin K.carrier ({y} : Set K.closure)) = σ := hτ
    rw [h3, h1] at h2
    simpa using h2.symm
  refine ⟨τ, fun z hz => ?_⟩
  obtain ⟨i, -, hi⟩ := hζ0.eq_pow_of_pow_eq_one hz
  rw [← hi, map_pow, hτζ, ← pow_mul, ← pow_mul, Nat.mul_comm]

/-- **惰性群は正規部分群**——`K^ur/K` が normal なので、`I_K` は制限射の核。

★`Found/PGC/InertiaIdentification.lean::normal_inertia` と同じ内容だが、
あちらは `inertia`(`Interface` 経由)で書かれていて循環するので、
ここでは `IntermediateField.restrictNormalHom_ker` から直接出す。

退化の自己検査:`unramifiedClosure K` を一般の中間体に替えると偽になる
(normal でない中間体の固定部分群は正規でない)。 -/
theorem normal_fixingSubgroup_unramifiedClosure (K : PAdicLocalField p) :
    ((unramifiedClosure K).fixingSubgroup).Normal := by
  haveI : Normal K.carrier (unramifiedClosure K) := instNormalUnramifiedClosure K
  rw [← IntermediateField.restrictNormalHom_ker (unramifiedClosure K)]
  exact MonoidHom.normal_ker _

/-- **Kummer 指標の乗法性**——`σ ↦ σβ/β` は `I_K` 上の準同型。

`σ₂β/β ∈ μ_n ⊆ K^ur`(在庫 `mem_unramifiedClosure_of_pow_eq_one`)なので
`σ₁` で固定される、というだけ。★`X^n − b` の既約性は使わない。

退化の自己検査:`hn` を落とすと `μ_n ⊆ K^ur` が崩れて `σ₁` が
`σ₂β/β` を固定しなくなるので偽。`hbE : β^n ∈ K^ur` を落とすと
`σβ/β` が `μ_n` に入らず、やはり偽。 -/
theorem kummer_mul (K : PAdicLocalField p) {n : ℕ} (hn : ¬ p ∣ n) {β : K.closure}
    (hβ0 : β ≠ 0) (hbE : β ^ n ∈ unramifiedClosure K)
    {σ₁ σ₂ : K.absGal} (h₁ : σ₁ ∈ (unramifiedClosure K).fixingSubgroup)
    (h₂ : σ₂ ∈ (unramifiedClosure K).fixingSubgroup) :
    (σ₁ * σ₂) β / β = (σ₁ β / β) * (σ₂ β / β) := by
  have hn0 : n ≠ 0 := by rintro rfl; exact hn (dvd_zero p)
  have hpow : ∀ σ : K.absGal, σ ∈ (unramifiedClosure K).fixingSubgroup → (σ β / β) ^ n = 1 := by
    intro σ hσ
    rw [div_pow, ← map_pow]
    rw [(IntermediateField.mem_fixingSubgroup_iff _ _).mp hσ _ hbE]
    exact div_self (pow_ne_zero n hβ0)
  have hfix : σ₁ (σ₂ β / β) = σ₂ β / β := by
    refine (IntermediateField.mem_fixingSubgroup_iff _ _).mp h₁ _ ?_
    exact mem_unramifiedClosure_of_pow_eq_one K (Nat.pos_of_ne_zero hn0) hn (hpow σ₂ h₂)
  have h : (σ₁ * σ₂) β = σ₁ (σ₂ β) := rfl
  rw [h]
  nth_rewrite 1 [show σ₂ β = (σ₂ β / β) * β from (div_mul_cancel₀ _ hβ0).symm]
  rw [map_mul, hfix]
  field_simp

/-- Kummer 指標の冪(`kummer_mul` の帰納)。 -/
theorem kummer_pow (K : PAdicLocalField p) {n : ℕ} (hn : ¬ p ∣ n) {β : K.closure}
    (hβ0 : β ≠ 0) (hbE : β ^ n ∈ unramifiedClosure K)
    {σ : K.absGal} (hσ : σ ∈ (unramifiedClosure K).fixingSubgroup) (m : ℕ) :
    (σ ^ m) β / β = (σ β / β) ^ m := by
  induction m with
  | zero => simpa using div_self hβ0
  | succ k ih =>
      rw [pow_succ, kummer_mul K hn hβ0 hbE (pow_mem hσ k) hσ, ih, pow_succ]

/-- **★★★★★★★★★★★★★★★★★★★★(F2)`Γ_K` の連続指標を惰性群に制限すると
`q − 1` 乗で消える**。

経路 C の (C-q) の上界の核心。証明の骨格はファイル冒頭の docstring を参照。

★要点は「`τβ/β ∈ K^ur`」——`‖τb/b‖ = 1` と在庫
`exists_pow_eq_mem_unramifiedClosure` から出るので、**`K^ur` の付値群を
調べる必要が無い**。当初の設計(`(K^ur)^× ≅ π^ℤ × 𝒪^×` を作る)より軽い。

退化の自己検査。

* `hn : ¬ p ∣ n` を落とすと**偽**。`K = ℚ_p`(`q = p`)、`n = p` で、
  `ℤ_p^× ≅ ℤ/(p−1) × ℤ_p` から `ℤ/p` への全射を局所類体論で `Γ_K` の指標に
  持ち上げると、それは分岐していて `f σ ≠ 1`(`σ ∈ I_K`)なのに
  `f σ ^ (p−1) = 1` は `f σ` の位数が `p` と `p−1` の両方を割ることを強いる。
* `hσ : σ ∈ I_K` を落とすと**偽**。次数 `n` の不分岐拡大が与える指標
  `Γ_K ↠ Gal(K_n/K) ≅ ℤ/n`(在庫 `exists_surjective_absGal_to_zmod`)は
  Frobenius 上で位数 `n` の値を取るので、`n ∤ q − 1` なら結論は成り立たない。
* `hf`(核が開)を落とすと本証明は届かない——Kummer 双対の全射性
  (`range_kummerMap`)が開核を使うため。主張自体の真偽は「`Γ_K` の
  有限指数部分群はすべて開か」(Nikolov–Segal 型の定理)に依存し、
  この木では判定していない。
* `n ∣ q − 1` のときは結論が**自明**になる(`ZMod n` 全体が `q−1` 乗で消える)。
  非自明なのは `n ∤ q − 1` の場合で、経路 C が使うのはまさにそこである。 -/
theorem pow_residueCard_sub_one_restrict (K : PAdicLocalField p) {n : ℕ} (hn : ¬ p ∣ n)
    (f : K.absGal →* Multiplicative (ZMod n))
    (hf : IsOpen ((MonoidHom.ker f : Subgroup K.absGal) : Set K.absGal))
    {σ : K.absGal} (hσ : σ ∈ (unramifiedClosure K).fixingSubgroup) :
    f σ ^ (Nat.card 𝓀[K.carrier] - 1) = 1 := by
  have hn0 : n ≠ 0 := by rintro rfl; exact hn (dvd_zero p)
  set q := Nat.card 𝓀[K.carrier] with hqdef
  have hq2 : 2 ≤ q := by
    have h1 : 1 < Fintype.card 𝓀[K.carrier] := Fintype.one_lt_card
    rw [hqdef, Nat.card_eq_fintype_card]; omega
  obtain ⟨β, hβ0, hbE, hkey⟩ := exists_kummer_root_inertia K hn f hf
  obtain ⟨τ, hτ⟩ := exists_absGal_pow_residueCard K hn
  set b := β ^ n with hbdef
  have hb0 : b ≠ 0 := pow_ne_zero n hβ0
  haveI : Normal K.carrier (unramifiedClosure K) := instNormalUnramifiedClosure K
  have hτbE : τ b ∈ unramifiedClosure K := by
    have h := AlgEquiv.restrictNormal_commutes τ (unramifiedClosure K) ⟨b, hbE⟩
    exact h ▸ (AlgEquiv.restrictNormal τ (unramifiedClosure K) ⟨b, hbE⟩).2
  have hτb0 : τ b ≠ 0 := by simpa using hb0
  have hdiv0 : τ b / b ≠ 0 := div_ne_zero hτb0 hb0
  have hdivE : τ b / b ∈ unramifiedClosure K := div_mem hτbE hbE
  have hnorm : ‖τ b / b‖ = 1 := by
    have hnb : ‖τ b‖ = ‖b‖ := norm_algEquiv_eq K τ b
    rw [norm_div, hnb, div_self (norm_ne_zero_iff.mpr hb0)]
  obtain ⟨w, hwE, hw⟩ := exists_pow_eq_mem_unramifiedClosure K hn hdivE hnorm
  have hw0 : w ≠ 0 := by
    intro h0; rw [h0, zero_pow hn0] at hw; exact hdiv0 hw.symm
  set v := τ β / β with hvdef
  have hvn : v ^ n = τ b / b := by rw [hvdef, div_pow, ← map_pow, hbdef]
  have hvw : (v / w) ^ n = 1 := by
    rw [div_pow, hvn, hw, div_self hdiv0]
  have hvwE : v / w ∈ unramifiedClosure K :=
    mem_unramifiedClosure_of_pow_eq_one K (Nat.pos_of_ne_zero hn0) hn hvw
  have hvE : v ∈ unramifiedClosure K := by
    rw [show v = (v / w) * w from (div_mul_cancel₀ _ hw0).symm]
    exact mul_mem hvwE hwE
  haveI := normal_fixingSubgroup_unramifiedClosure K
  have hσ'mem : τ⁻¹ * σ * τ ∈ (unramifiedClosure K).fixingSubgroup := by
    have h := (normal_fixingSubgroup_unramifiedClosure K).conj_mem σ hσ τ⁻¹
    simpa using h
  have hkap : σ (τ β) = (σ β / β) * (τ β) := by
    have hfixv : σ v = v := (IntermediateField.mem_fixingSubgroup_iff _ _).mp hσ _ hvE
    have h1 : τ β = v * β := by rw [hvdef]; field_simp
    rw [h1, map_mul, hfixv]
    field_simp
  have hconj : (τ⁻¹ * σ * τ) β / β = τ⁻¹ (σ β / β) := by
    have h2 : (τ⁻¹ * σ * τ) β = τ⁻¹ (σ (τ β)) := rfl
    rw [h2, hkap, map_mul]
    have h3 : τ⁻¹ (τ β) = β := by simp
    rw [h3]
    field_simp
  set κ := σ β / β with hκdef
  have hκn : κ ^ n = 1 := by
    rw [hκdef, div_pow, ← map_pow, ← hbdef,
      (IntermediateField.mem_fixingSubgroup_iff _ _).mp hσ _ hbE, div_self hb0]
  have hκ'n : (τ⁻¹ κ) ^ n = 1 := by rw [← map_pow, hκn, map_one]
  have hτκ' : τ (τ⁻¹ κ) = (τ⁻¹ κ) ^ q := hτ _ hκ'n
  have hκeq : κ = (τ⁻¹ κ) ^ q := by rw [← hτκ']; simp
  have hkey2 : ((τ⁻¹ * σ * τ) ^ q) β / β = σ β / β := by
    rw [kummer_pow K hn hβ0 hbE hσ'mem q, hconj, ← hκdef, ← hκeq]
  have hfeq : f ((τ⁻¹ * σ * τ) ^ q) = f σ := hkey _ _ (pow_mem hσ'mem q) hσ hkey2
  have hfσ' : f (τ⁻¹ * σ * τ) = f σ := by
    rw [map_mul, map_mul, map_inv, mul_right_comm, inv_mul_cancel, one_mul]
  rw [map_pow, hfσ'] at hfeq
  have h5 : f σ ^ (q - 1) * f σ = 1 * f σ := by
    rw [one_mul, ← pow_succ, show q - 1 + 1 = q by omega]
    exact hfeq
  exact mul_right_cancel h5


/-! ## `K^ur` の付値群は `K` の付値群 —— 不分岐なら `e = 1`

(F1) `Hom_cont(I_K, ℤ/n) ≅ ℤ/n` を出すには、`(K^ur)^×` を `n` 乗で割った群が
**ちょうど位数 `n`** であることが要る。上界(F2)では避けられた
「`K^ur` の付値群」を、ここで初めて調べる。

要点は 3 行で済む:`𝒪_{K(y)}` は離散付値環(在庫)なので素元 `ϖ` を持ち、
`K` の素元 `π` は `π = u·ϖ^j` と書ける。`Ideal.span {π} ≤ 𝔪^j` から
`j ≤ e = 1`(不分岐)、`‖π‖ < 1` から `j ≥ 1`。よって `‖π‖ = ‖ϖ‖`。
-/

/-- **不分岐拡大では `K` の素元が `𝒪_{K(y)}` の素元**(ノルムで見た形)。

退化の自己検査:`hu : IsUnramifiedAdjoin K y` を落とすと**偽**。完全分岐な
`K(y) = K(π^{1/e})` では `‖ϖ‖ = ‖π‖^{1/e} ≠ ‖π‖`。`hπ` / `hϖ` の既約性を
落とすと `j` の決定ができず、主張自体が意味を失う。 -/
theorem norm_irreducible_adjoinIntegers (K : PAdicLocalField p) {y : K.closure}
    (hu : IsUnramifiedAdjoin K y) {π : 𝒪[K.carrier]} (hπ : Irreducible π)
    {ϖ : adjoinIntegers K y} (hϖ : Irreducible ϖ) :
    ‖(ϖ : IntermediateField.adjoin K.carrier ({y} : Set K.closure))‖ = ‖(π : K.carrier)‖ := by
  haveI := isDiscreteValuationRing_adjoinIntegers K y
  haveI := isDiscreteValuationRing_carrierIntegers K
  have hπlt : ‖(π : K.carrier)‖ < 1 := by
    have hle : ‖(π : K.carrier)‖ ≤ 1 := by
      have h := π.2; rw [Valued.integer.mem_iff] at h; exact h
    have hne : ‖(π : K.carrier)‖ ≠ 1 := fun h =>
      hπ.not_isUnit (Valued.integer.isUnit_iff_norm_eq_one.mpr h)
    exact lt_of_le_of_ne hle hne
  have hπ0 : ‖(π : K.carrier)‖ ≠ 0 := by
    intro h
    have h2 : (π : K.carrier) = 0 := by simpa using h
    exact hπ.ne_zero (Subtype.ext h2)
  set ip := algebraMap 𝒪[K.carrier] (adjoinIntegers K y) π with hip
  have hipn : ‖((ip : adjoinIntegers K y) :
      IntermediateField.adjoin K.carrier ({y} : Set K.closure))‖ = ‖(π : K.carrier)‖ :=
    norm_algebraMap_adjoin K y _
  have hip0 : ip ≠ 0 := by
    intro h0
    rw [h0] at hipn
    simp only [ZeroMemClass.coe_zero, norm_zero] at hipn
    exact hπ0 hipn.symm
  obtain ⟨j, u, hju⟩ := IsDiscreteValuationRing.eq_unit_mul_pow_irreducible hip0 hϖ
  have hun : ‖((u : adjoinIntegers K y) :
      IntermediateField.adjoin K.carrier ({y} : Set K.closure))‖ = 1 :=
    (isUnit_adjoinIntegers_iff K y _).mp u.isUnit
  have hnormj : ‖(π : K.carrier)‖
      = ‖(ϖ : IntermediateField.adjoin K.carrier ({y} : Set K.closure))‖ ^ j := by
    rw [← hipn, hju]
    push_cast
    rw [norm_mul, hun, one_mul, norm_pow]
  have hjle : j ≤ 1 := by
    have hmap : Ideal.map (algebraMap 𝒪[K.carrier] (adjoinIntegers K y))
        (IsLocalRing.maximalIdeal 𝒪[K.carrier])
        ≤ (IsLocalRing.maximalIdeal (adjoinIntegers K y)) ^ j := by
      rw [(IsDiscreteValuationRing.irreducible_iff_uniformizer π).mp hπ,
        (IsDiscreteValuationRing.irreducible_iff_uniformizer ϖ).mp hϖ,
        Ideal.map_span, Set.image_singleton, Ideal.span_singleton_pow,
        Ideal.span_singleton_le_span_singleton]
      exact ⟨u, by rw [← hip, hju]; ring⟩
    have h := le_ramificationIdx_of_map_le_pow (ramificationIndex_ne_zero K y) hmap
    rw [show Ideal.ramificationIdx (IsLocalRing.maximalIdeal 𝒪[K.carrier])
        (IsLocalRing.maximalIdeal (adjoinIntegers K y)) = ramificationIndex K y from rfl,
      hu] at h
    exact h
  have hjge : 1 ≤ j := by
    rcases Nat.eq_zero_or_pos j with h0 | h1
    · rw [h0, pow_zero] at hnormj; exact absurd hnormj (by linarith)
    · exact h1
  have hj : j = 1 := le_antisymm hjle hjge
  rw [hj, pow_one] at hnormj
  exact hnormj.symm

/-- `K^ur` の整な元のノルムは `‖π‖^k`(`k : ℕ`)。 -/
theorem exists_norm_eq_pow_of_le_one (K : PAdicLocalField p) {π : 𝒪[K.carrier]}
    (hπ : Irreducible π) {x : K.closure} (hx : x ∈ unramifiedClosure K) (hx0 : x ≠ 0)
    (hle : ‖x‖ ≤ 1) : ∃ k : ℕ, ‖x‖ = ‖(π : K.carrier)‖ ^ k := by
  obtain ⟨y, huy, hxy⟩ := (mem_unramifiedClosure_iff K x).mp hx
  haveI := isDiscreteValuationRing_adjoinIntegers K y
  have hxL : ‖(⟨x, hxy⟩ : IntermediateField.adjoin K.carrier ({y} : Set K.closure))‖ = ‖x‖ := rfl
  set z : adjoinIntegers K y := ⟨⟨x, hxy⟩,
    (by rw [hxL]; exact hle :
      ‖(⟨x, hxy⟩ : IntermediateField.adjoin K.carrier ({y} : Set K.closure))‖ ≤ 1)⟩ with hz
  have hz0 : z ≠ 0 := by
    intro h
    apply hx0
    have hzx : ((z : IntermediateField.adjoin K.carrier ({y} : Set K.closure)) : K.closure)
      = x := rfl
    rw [← hzx, h]
    rfl
  obtain ⟨ϖ, hϖ⟩ := IsDiscreteValuationRing.exists_irreducible (adjoinIntegers K y)
  obtain ⟨k, u, hku⟩ := IsDiscreteValuationRing.eq_unit_mul_pow_irreducible hz0 hϖ
  have hun : ‖((u : adjoinIntegers K y) :
      IntermediateField.adjoin K.carrier ({y} : Set K.closure))‖ = 1 :=
    (isUnit_adjoinIntegers_iff K y _).mp u.isUnit
  refine ⟨k, ?_⟩
  have h1 : ‖(z : IntermediateField.adjoin K.carrier ({y} : Set K.closure))‖ = ‖x‖ := hxL
  rw [← h1, hku]
  push_cast
  rw [norm_mul, hun, one_mul, norm_pow, norm_irreducible_adjoinIntegers K huy hπ hϖ]

/-- **★★★★★★★★`K^ur` の付値群は `K` の付値群**——`‖x‖ = ‖π‖^k`(`k : ℤ`)。

退化の自己検査:`hx : x ∈ K^ur` を落とすと**偽**(`x = π^{1/2}` のノルムは
`‖π‖^{1/2}` で `‖π‖^ℤ` に無い)。`hx0` を落とすと `‖0‖ = 0` はどの `‖π‖^k` とも
等しくないので偽。 -/
theorem exists_norm_eq_zpow (K : PAdicLocalField p) {π : 𝒪[K.carrier]}
    (hπ : Irreducible π) {x : K.closure} (hx : x ∈ unramifiedClosure K) (hx0 : x ≠ 0) :
    ∃ k : ℤ, ‖x‖ = ‖(π : K.carrier)‖ ^ k := by
  rcases le_or_gt ‖x‖ 1 with hle | hgt
  · obtain ⟨k, hk⟩ := exists_norm_eq_pow_of_le_one K hπ hx hx0 hle
    exact ⟨(k : ℤ), by rw [hk, zpow_natCast]⟩
  · have hxi : x⁻¹ ∈ unramifiedClosure K := inv_mem hx
    have hxi0 : x⁻¹ ≠ 0 := inv_ne_zero hx0
    have hxile : ‖x⁻¹‖ ≤ 1 := by
      rw [norm_inv]
      exact inv_le_one_of_one_le₀ hgt.le
    obtain ⟨k, hk⟩ := exists_norm_eq_pow_of_le_one K hπ hxi hxi0 hxile
    refine ⟨-(k : ℤ), ?_⟩
    rw [zpow_neg, zpow_natCast, ← hk, norm_inv, inv_inv]

/-! ## `(K^ur)^× / ((K^ur)^×)^n` は位数 `n` の巡回群 -/

/-- **★★★★★★★★`(K^ur)^×` を `n` 乗で割ると `[π]` の生成する位数 `n` の巡回群**。

* 生成:`‖a‖ = ‖π‖^k` なので `a·π^{-k}` はノルム `1`、よって在庫
  `exists_pow_eq_mem_unramifiedClosure` で `n` 乗数。
* 位数:`π^m = b^n` なら `‖π‖^m = ‖π‖^{jn}`、`0 < ‖π‖ < 1` から `m = jn`。

退化の自己検査:`hn : ¬ p ∣ n` を落とすと**偽**。`n = p` では
`𝒪_{K^ur}^×` が `p` 可除でないので商は `ℤ/p` より真に大きい。 -/
theorem exists_generator_quotient_powRange (K : PAdicLocalField p) {n : ℕ} (hn : ¬ p ∣ n) :
    ∃ g : ((↥(unramifiedClosure K))ˣ ⧸
        (powMonoidHom n : (↥(unramifiedClosure K))ˣ →* (↥(unramifiedClosure K))ˣ).range),
      (∀ q, q ∈ Subgroup.zpowers g) ∧ orderOf g = n := by
  have hn0 : n ≠ 0 := by rintro rfl; exact hn (dvd_zero p)
  haveI := isDiscreteValuationRing_carrierIntegers K
  set E := unramifiedClosure K with hE
  set H := (powMonoidHom n : (↥E)ˣ →* (↥E)ˣ).range with hH
  obtain ⟨π, hπ⟩ := IsDiscreteValuationRing.exists_irreducible 𝒪[K.carrier]
  have hπlt : ‖(π : K.carrier)‖ < 1 := by
    have hle : ‖(π : K.carrier)‖ ≤ 1 := by
      have h := π.2; rw [Valued.integer.mem_iff] at h; exact h
    have hne : ‖(π : K.carrier)‖ ≠ 1 := fun h =>
      hπ.not_isUnit (Valued.integer.isUnit_iff_norm_eq_one.mpr h)
    exact lt_of_le_of_ne hle hne
  have hπ0 : (π : K.carrier) ≠ 0 := fun h => hπ.ne_zero (Subtype.ext h)
  have hπpos : 0 < ‖(π : K.carrier)‖ := norm_pos_iff.mpr hπ0
  have hcoez : ∀ (a : (↥E)ˣ) (m : ℤ),
      (((a ^ m : (↥E)ˣ) : ↥E) : K.closure) = ((a : ↥E) : K.closure) ^ m := by
    intro a m
    rw [Units.val_zpow_eq_zpow_val]
    exact map_zpow₀ E.val _ m
  have hne0 : ∀ a : (↥E)ˣ, ((a : ↥E) : K.closure) ≠ 0 := by
    intro a h
    exact a.ne_zero (Subtype.ext h)
  have hval : ∀ a : (↥E)ˣ, ∃ k : ℤ,
      ‖((a : ↥E) : K.closure)‖ = ‖(π : K.carrier)‖ ^ k := fun a =>
    exists_norm_eq_zpow K hπ (a : ↥E).2 (hne0 a)
  set πc : K.closure := algebraMap K.carrier K.closure (π : K.carrier) with hπc
  have hπcnorm : ‖πc‖ = ‖(π : K.carrier)‖ := spectralNorm_extends _
  have hπc0 : πc ≠ 0 := by
    rw [hπc, Ne, map_eq_zero_iff _ (algebraMap K.carrier K.closure).injective]
    exact hπ0
  set πF : (↥E)ˣ := Units.mk0 (⟨πc, E.algebraMap_mem _⟩ : ↥E)
    (fun h => hπc0 (congrArg Subtype.val h)) with hπF
  have hπFc : ((πF : ↥E) : K.closure) = πc := rfl
  refine ⟨QuotientGroup.mk' H πF, ?_, ?_⟩
  · intro q
    induction q using QuotientGroup.induction_on with
    | _ a =>
      obtain ⟨k, hk⟩ := hval a
      set c : (↥E)ˣ := a * πF ^ (-k) with hc
      have hcnorm : ‖((c : ↥E) : K.closure)‖ = 1 := by
        have h1 : ((c : ↥E) : K.closure)
            = ((a : ↥E) : K.closure) * ((πF : ↥E) : K.closure) ^ (-k) := by
          rw [← hcoez πF (-k), hc]
          push_cast
          ring
        rw [h1, norm_mul, norm_zpow, hk, hπFc, hπcnorm, ← zpow_add₀ (ne_of_gt hπpos)]
        simp
      obtain ⟨w, hwE, hw⟩ := exists_pow_eq_mem_unramifiedClosure K hn (c : ↥E).2 hcnorm
      have hw0 : w ≠ 0 := by
        intro h0
        rw [h0, zero_pow hn0] at hw
        exact hne0 c hw.symm
      set wU : (↥E)ˣ := Units.mk0 (⟨w, hwE⟩ : ↥E) (fun h => hw0 (congrArg Subtype.val h)) with hwU
      have hcH : c ∈ H := by
        refine ⟨wU, ?_⟩
        apply Units.ext
        apply Subtype.ext
        show w ^ n = ((c : ↥E) : K.closure)
        exact hw
      have hca : a = c * πF ^ k := by rw [hc]; group
      refine ⟨k, ?_⟩
      show (QuotientGroup.mk' H πF) ^ k = QuotientGroup.mk' H a
      rw [(map_zpow (QuotientGroup.mk' H) πF k).symm]
      show ((πF ^ k : (↥E)ˣ) : (↥E)ˣ ⧸ H) = (a : (↥E)ˣ ⧸ H)
      rw [QuotientGroup.eq]
      have h9 : (πF ^ k)⁻¹ * a = c := by
        rw [hca, mul_comm c (πF ^ k), inv_mul_cancel_left]
      rw [h9]
      exact hcH
  · have hgn : (QuotientGroup.mk' H πF) ^ n = 1 := by
      rw [← map_pow]
      exact (QuotientGroup.eq_one_iff _).mpr ⟨πF, rfl⟩
    have hdvd1 : orderOf (QuotientGroup.mk' H πF) ∣ n := orderOf_dvd_of_pow_eq_one hgn
    have hdvd2 : n ∣ orderOf (QuotientGroup.mk' H πF) := by
      set m := orderOf (QuotientGroup.mk' H πF) with hm
      have hpm : (QuotientGroup.mk' H πF) ^ m = 1 := pow_orderOf_eq_one _
      have h1 : πF ^ m ∈ H := by
        rw [← QuotientGroup.eq_one_iff, ← QuotientGroup.mk'_apply, map_pow]
        exact hpm
      obtain ⟨b, hb⟩ := h1
      obtain ⟨j, hj⟩ := hval b
      have hbn : ((b : ↥E) : K.closure) ^ n = ((πF : ↥E) : K.closure) ^ m := by
        have hcong := congrArg (fun t : (↥E)ˣ => ((t : ↥E) : K.closure)) hb
        simpa [powMonoidHom] using hcong
      have h2 : ‖(π : K.carrier)‖ ^ (j * n) = ‖(π : K.carrier)‖ ^ (m : ℤ) := by
        have hl : ‖((b : ↥E) : K.closure) ^ n‖ = ‖(π : K.carrier)‖ ^ (j * n) := by
          rw [norm_pow, hj, ← zpow_natCast (‖(π : K.carrier)‖ ^ j) n, ← zpow_mul]
        have hr : ‖((πF : ↥E) : K.closure) ^ m‖ = ‖(π : K.carrier)‖ ^ (m : ℤ) := by
          rw [norm_pow, hπFc, hπcnorm, zpow_natCast]
        rw [← hl, ← hr, hbn]
      have h3 : j * n = (m : ℤ) :=
        zpow_right_injective₀ hπpos (ne_of_lt hπlt) h2
      exact Int.natCast_dvd_natCast.mp ⟨j, by linarith⟩
    exact Nat.dvd_antisymm hdvd1 hdvd2

/-! ## (F1)`Hom_cont(I_K, ℤ/n)` は位数 `n` の巡回群 -/

/-- 位相群の同型に沿った `contHom` の**乗法的**同型
(`Found/PGC/ContinuousHomCount.lean::contHomEquiv` は `Equiv` 止まり)。 -/
noncomputable def contHomMulEquiv {G H : Type*} [Group G] [TopologicalSpace G]
    [SeparatelyContinuousMul G] [Group H] [TopologicalSpace H] [SeparatelyContinuousMul H]
    (α : ContinuousMulEquiv G H) (A : Type*) [CommGroup A] :
    contHom H A ≃* contHom G A where
  toEquiv := contHomEquiv α A
  map_mul' _ _ := rfl

open KummerDual in
/-- 係数群を取り替える `contHom` の**乗法的**同型。 -/
noncomputable def contHomMulEquivOfMulEquiv {G : Type*} [Group G] [TopologicalSpace G]
    [SeparatelyContinuousMul G] {A B : Type*} [CommGroup A] [CommGroup B] (e : A ≃* B) :
    contHom G A ≃* contHom G B where
  toEquiv := contHomEquivOfMulEquiv e
  map_mul' f g := by
    apply Subtype.ext
    refine MonoidHom.ext fun x => ?_
    show e ((f : G →* A) x * (g : G →* A) x) = e ((f : G →* A) x) * e ((g : G →* A) x)
    exact map_mul e _ _

open KummerDual in
/-- **`Hom_cont(I_K, ℤ/n) ≅ (K^ur)^× / ((K^ur)^×)^n`**——`K^ur` 上の Kummer 双対を
`I_K` に載せ替えたもの。5 段の合成:

1. `fixingSubgroupContinuousMulEquivInf`(第 1022、無限次でも通る)
2. `autCongrContinuousMulEquiv`(第 993、代数閉包の取り替え)
3. 係数群 `ℤ/n ≃* μ_n(Ω)`(`zmodMulEquivRootsOfUnity`)
4. `range_kummerMap`(第 1020、Kummer 双対の本体)
5. 第一同型定理と `ker_kummerMap` -/
theorem nonempty_contHom_inertia_mulEquiv (K : PAdicLocalField p) {n : ℕ} (hn : ¬ p ∣ n) :
    Nonempty (contHom ((unramifiedClosure K).fixingSubgroup) (Multiplicative (ZMod n)) ≃*
      ((↥(unramifiedClosure K))ˣ ⧸
        (powMonoidHom n : (↥(unramifiedClosure K))ˣ →* (↥(unramifiedClosure K))ˣ).range)) := by
  have hn0 : n ≠ 0 := by rintro rfl; exact hn (dvd_zero p)
  set E := unramifiedClosure K with hE
  obtain ⟨ζ0, hmem, hζ0⟩ := exists_isPrimitiveRoot_mem_unramifiedClosure K n
    (Nat.pos_of_ne_zero hn0) hn
  have hζ : IsPrimitiveRoot (⟨ζ0, hmem⟩ : ↥E) n :=
    IsPrimitiveRoot.of_map_of_injective (f := algebraMap ↥E K.closure) hζ0
      (algebraMap ↥E K.closure).injective
  refine ⟨?_⟩
  refine MulEquiv.trans ?_ (QuotientGroup.quotientMulEquivOfEq (ker_kummerMap hζ hn0))
  refine MulEquiv.trans (contHomMulEquiv (fixingSubgroupContinuousMulEquivInf E)
    (Multiplicative (ZMod n))).symm ?_
  refine MulEquiv.trans (contHomMulEquiv
    (autCongrContinuousMulEquiv (IsAlgClosure.equiv ↥E (AlgebraicClosure ↥E) K.closure))
    (Multiplicative (ZMod n))) ?_
  refine MulEquiv.trans (contHomMulEquivOfMulEquiv (zmodMulEquivRootsOfUnity hζ hn0)) ?_
  exact ((MulEquiv.subgroupCongr (range_kummerMap hζ hn0)).symm).trans
    (QuotientGroup.quotientKerEquivRange (kummerMap hζ hn0)).symm

/-- **★★★★★★★★★★★★★★★★(F1)`Hom_cont(I_K, ℤ/n)` は位数 `n` の巡回群**
(`p ∤ n`)。

`K^ur` 上の Kummer 理論:`Hom_cont(I_K, μ_n) ≅ (K^ur)^×/n ≅ ℤ/n`。
右の同型は「`𝒪_{K^ur}^×` は `n` 可除(在庫)」と「`K^ur` の付値群は `K` の付値群」
(`exists_norm_eq_zpow`)から。

退化の自己検査。

* `hn : ¬ p ∣ n` を落とすと**偽**。`n = p` では `Hom_cont(I_K, ℤ/p)` は
  野性的惰性群のぶんだけ大きく、位数 `p` にならない。
* `I_K` を `Γ_K` に替えると**偽**——`N_n(Γ_K) = n·gcd(n,q−1)` は一般に `n` より大きい
  (在庫 `contHomCard_absGal_of_dvd`:`n ∣ q−1` なら `n²`)。
* `n = 1` では両辺が自明群になるが、主張は真のまま(退化して自明)。 -/
theorem contHom_fixingSubgroup_unramifiedClosure (K : PAdicLocalField p) {n : ℕ} (hn : ¬ p ∣ n) :
    Nonempty (contHom ((unramifiedClosure K).fixingSubgroup) (Multiplicative (ZMod n))
      ≃* Multiplicative (ZMod n)) := by
  obtain ⟨φ⟩ := nonempty_contHom_inertia_mulEquiv K hn
  obtain ⟨g, hgen, hord⟩ := exists_generator_quotient_powRange K hn
  haveI hc : IsCyclic ((↥(unramifiedClosure K))ˣ ⧸
      (powMonoidHom n : (↥(unramifiedClosure K))ˣ →* (↥(unramifiedClosure K))ˣ).range) :=
    ⟨⟨g, hgen⟩⟩
  have hcard : Nat.card ((↥(unramifiedClosure K))ˣ ⧸
      (powMonoidHom n : (↥(unramifiedClosure K))ˣ →* (↥(unramifiedClosure K))ˣ).range) = n := by
    have h1 : Nat.card ((↥(unramifiedClosure K))ˣ ⧸
        (powMonoidHom n : (↥(unramifiedClosure K))ˣ →* (↥(unramifiedClosure K))ˣ).range)
        = orderOf g := by
      rw [← Nat.card_zpowers g]
      exact (Nat.card_congr (Equiv.subtypeUnivEquiv hgen)).symm
    rw [h1, hord]
  refine ⟨φ.trans ?_⟩
  have h2 := (zmodCyclicMulEquiv hc).symm
  rw [hcard] at h2
  exact h2


/-! ## (F1)+(F2) の消費 —— 惰性群への制限射の像は `gcd(n, q−1)` 個以下

(C-q) の上界 `N_n(Γ_K) ≤ n · gcd(n, q−1)`(ノード G2)は、連続指標の群を
惰性群への制限射

  `restrictInertia : Hom_cont(Γ_K, ℤ/n) → Hom_cont(I_K, ℤ/n)`

で切って `#核 × #像` に分けることで出る。本節はそのうち**像の側**を
(F1)(F2)から完成させる。★核の側(`#核 ≤ n`、すなわち F3
`Hom_cont(Gal(K^ur/K), ℤ/n)` の個数)は本ファイルでは埋まっていない
——`contHomCard_absGal_le_of_card_ker` がその一点だけを仮定に残した形で
G2 を用意してある。
-/

/-- 惰性群への制限(連続指標の群としての準同型)。 -/
noncomputable def restrictInertia (K : PAdicLocalField p) (n : ℕ) :
    contHom K.absGal (Multiplicative (ZMod n)) →*
      contHom ((unramifiedClosure K).fixingSubgroup) (Multiplicative (ZMod n)) where
  toFun f := ⟨(f : K.absGal →* Multiplicative (ZMod n)).comp
      ((unramifiedClosure K).fixingSubgroup).subtype,
    comp_mem_contHom _ continuous_subtype_val f.2⟩
  map_one' := rfl
  map_mul' _ _ := rfl

open KummerDual in
/-- **`Hom_cont(I_K, ℤ/n)` の `m` 乗根の個数はちょうど `gcd(n, m)`**。

(F1)で `Hom_cont(I_K, ℤ/n)` が位数 `n` の巡回群だと判っているので、
`KummerDuality` の群論補題 `card_ker_powMonoidHom` がそのまま当たる。

退化の自己検査:`hn : ¬ p ∣ n` を落とすと(F1)が崩れるので**偽**。
`m = 0` では `gcd(n,0) = n` で「全体」を表し、主張は真のまま(自明側へ退化)。 -/
theorem card_ker_powMonoidHom_contHom_inertia (K : PAdicLocalField p) {n : ℕ} (hn : ¬ p ∣ n)
    (m : ℕ) :
    Nat.card ((powMonoidHom m :
        contHom ((unramifiedClosure K).fixingSubgroup) (Multiplicative (ZMod n)) →*
        contHom ((unramifiedClosure K).fixingSubgroup) (Multiplicative (ZMod n))).ker)
      = Nat.gcd n m := by
  have hn0 : n ≠ 0 := by rintro rfl; exact hn (dvd_zero p)
  haveI : NeZero n := ⟨hn0⟩
  obtain ⟨φ⟩ := contHom_fixingSubgroup_unramifiedClosure K hn
  haveI : Finite (contHom ((unramifiedClosure K).fixingSubgroup) (Multiplicative (ZMod n))) :=
    Finite.of_equiv _ φ.toEquiv.symm
  haveI : IsCyclic (contHom ((unramifiedClosure K).fixingSubgroup) (Multiplicative (ZMod n))) :=
    (MulEquiv.isCyclic φ).mpr inferInstance
  have hcard : Nat.card
      (contHom ((unramifiedClosure K).fixingSubgroup) (Multiplicative (ZMod n))) = n := by
    rw [Nat.card_congr φ.toEquiv]
    simp
  have hker : (powMonoidHom m :
      contHom ((unramifiedClosure K).fixingSubgroup) (Multiplicative (ZMod n)) →*
      contHom ((unramifiedClosure K).fixingSubgroup) (Multiplicative (ZMod n))).ker
      = (powMonoidHom (Nat.gcd n m) :
      contHom ((unramifiedClosure K).fixingSubgroup) (Multiplicative (ZMod n)) →*
      contHom ((unramifiedClosure K).fixingSubgroup) (Multiplicative (ZMod n))).ker := by
    ext x
    simp only [MonoidHom.mem_ker, powMonoidHom_apply]
    constructor
    · intro hx
      have h1 : orderOf x ∣ m := orderOf_dvd_of_pow_eq_one hx
      have h2 : orderOf x ∣ n := by
        have hd := orderOf_dvd_natCard x
        rwa [hcard] at hd
      exact orderOf_dvd_iff_pow_eq_one.mp (Nat.dvd_gcd h2 h1)
    · intro hx
      have h1 : orderOf x ∣ Nat.gcd n m := orderOf_dvd_of_pow_eq_one hx
      exact orderOf_dvd_iff_pow_eq_one.mp (h1.trans (Nat.gcd_dvd_right n m))
  rw [hker]
  exact card_ker_powMonoidHom (by rw [hcard]; exact Nat.gcd_dvd_left n m)

/-- **(F2)の言い換え**——制限射の像は `q − 1` 乗で消える部分群に入る。 -/
theorem range_restrictInertia_le (K : PAdicLocalField p) {n : ℕ} (hn : ¬ p ∣ n) :
    (restrictInertia K n).range ≤ (powMonoidHom (Nat.card 𝓀[K.carrier] - 1) :
      contHom ((unramifiedClosure K).fixingSubgroup) (Multiplicative (ZMod n)) →*
      contHom ((unramifiedClosure K).fixingSubgroup) (Multiplicative (ZMod n))).ker := by
  rintro _ ⟨f, rfl⟩
  rw [MonoidHom.mem_ker, powMonoidHom_apply]
  apply Subtype.ext
  refine MonoidHom.ext fun σ => ?_
  show (f : K.absGal →* Multiplicative (ZMod n)) (σ : K.absGal) ^
    (Nat.card 𝓀[K.carrier] - 1) = 1
  exact pow_residueCard_sub_one_restrict K hn _ f.2 σ.2

/-- **★★★★★★★★★★★★(F1)+(F2):制限射の像は `gcd(n, q−1)` 個以下**。

(C-q) の上界 `N_n(Γ_K) ≤ n · gcd(n, q−1)` の**像の側**。

退化の自己検査:`hn : ¬ p ∣ n` を落とすと(F1)(F2)がともに崩れて**偽**。
`n ∣ q − 1` のときは `gcd = n` で上界が `n` になり、下界(在庫 G1 の
`contHomCard_absGal_of_dvd` が使う `N_n = n²`)と整合する
——このとき像は `Hom_cont(I_K, ℤ/n)` 全体である。 -/
theorem card_range_restrictInertia_le (K : PAdicLocalField p) {n : ℕ} (hn : ¬ p ∣ n) :
    Nat.card ((restrictInertia K n).range)
      ≤ Nat.gcd n (Nat.card 𝓀[K.carrier] - 1) := by
  have hn0 : n ≠ 0 := by rintro rfl; exact hn (dvd_zero p)
  haveI : NeZero n := ⟨hn0⟩
  obtain ⟨φ⟩ := contHom_fixingSubgroup_unramifiedClosure K hn
  haveI : Finite (contHom ((unramifiedClosure K).fixingSubgroup) (Multiplicative (ZMod n))) :=
    Finite.of_equiv _ φ.toEquiv.symm
  have hle := range_restrictInertia_le K hn
  have h1 : Nat.card ((restrictInertia K n).range)
      ≤ Nat.card ((powMonoidHom (Nat.card 𝓀[K.carrier] - 1) :
        contHom ((unramifiedClosure K).fixingSubgroup) (Multiplicative (ZMod n)) →*
        contHom ((unramifiedClosure K).fixingSubgroup) (Multiplicative (ZMod n))).ker) :=
    Nat.card_le_card_of_injective (Subgroup.inclusion hle) (Subgroup.inclusion_injective hle)
  rw [card_ker_powMonoidHom_contHom_inertia K hn] at h1
  exact h1

/-- **(G2)の骨格**——`#ker(restrictInertia) ≤ n`(= F3、本ファイルでは未完)を
仮定すれば、`N_n(Γ_K) ≤ n · gcd(n, q−1)` が出る。

`#Hom_cont(Γ_K, ℤ/n) = #核 × #像` を第一同型定理で分け、像は
`card_range_restrictInertia_le`(本ファイルで完成)で押さえる。

★残っている一点は「`Hom_cont(Γ_K, ℤ/n)` のうち `I_K` 上自明なものは `n` 個以下」
——`Γ_K/I_K ≅ Gal(K^ur/K) ≅ Ẑ` の連続指標の個数であり、`Ẑ` の位相的生成元
(Frobenius)の存在、すなわち不分岐拡大の塔の射影極限が要る。
`ResearchPaper/pgc-goal.md` の経路 C ではノード F3 に当たる。

退化の自己検査:`hn` を落とすと像の評価(F1/F2)が崩れて**偽**。
`hker` を落とすと、核の側が抑えられないので上界は出ない
(`Γ_K` の連続指標は `I_K` 上自明なものだけで無限個ありうる、という
可能性を排除しているのが `hker` である)。 -/
theorem contHomCard_absGal_le_of_card_ker (K : PAdicLocalField p) {n : ℕ} (hn : ¬ p ∣ n)
    (hker : Nat.card ((restrictInertia K n).ker) ≤ n) :
    contHomCard K.absGal n ≤ n * Nat.gcd n (Nat.card 𝓀[K.carrier] - 1) := by
  have h1 : Nat.card ((restrictInertia K n).ker) * ((restrictInertia K n).ker).index
      = Nat.card (contHom K.absGal (Multiplicative (ZMod n))) := Subgroup.card_mul_index _
  have h2 : ((restrictInertia K n).ker).index = Nat.card ((restrictInertia K n).range) :=
    Nat.card_congr (QuotientGroup.quotientKerEquivRange _).toEquiv
  have h3 : contHomCard K.absGal n
      = Nat.card ((restrictInertia K n).ker) * Nat.card ((restrictInertia K n).range) := by
    rw [← h2, h1]
    rfl
  rw [h3]
  exact Nat.mul_le_mul hker (card_range_restrictInertia_le K hn)

end ABC3.Found.PGC
