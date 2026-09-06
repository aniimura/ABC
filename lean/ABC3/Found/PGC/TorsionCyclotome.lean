import ABC3.Found.PGC.TopAbelianization
import ABC3.Found.PGC.InertiaTransport

/-!
# 円分子 `Λ_n(H) = tors_{p^n}(H^{ab})` —— 定義・共役作用・移送

[pGC] Proposition 1.1 への**経路 Λ**の節点 Λ2。Λ1
(`Found/PGC/TopAbelianization.lean`)の上に乗り、Λ10
(`Found/PGC/CyclotomeTransport.lean`)がこれを消費する。

## 定義の設計——`H` を**データとして持たない**

素朴には「`Γ_K` の開正規部分群 `H` に対して `Λ_n(H)` を定義する」と書きたくなるが、
それは Proposition 1.2 の `∀ RD` と同じ罠(自由なパラメータが結論に現れる)を招く
(`Check/PGC/Prop12ForallRD.lean` の 10 例目)。本ファイルは

  `cyclotome A m := ker (x ↦ x^m : A^{ab} → A^{ab})`

と、**位相群 `A` だけ**を引数に取って定義する。`H` は `cyclotome ↥H (p^n)` のように
「`A := ↥H`」として現れるだけで、定義のデータではない。したがって
`H` は下流の結論の statement に決して現れない——現れるとすれば
`∃ H` / `∀ H` の内側だけである。

## 本ファイルが確立したもの

1. **`cyclotome A m`**(`A^{ab}` の `m` 捩れ)と `mem_cyclotome`。
2. **移送** `cyclotomeEquiv : (A ≃ₜ* B) → cyclotome A m ≃* cyclotome B m`、
   および関手性(`cyclotomeEquiv_trans` / `_refl` / `_congr`)。
3. **共役作用** `conjSubgroupCME g : ↥H ≃ₜ* ↥H`(`H` 正規)と、それが誘導する
   `cyclotomeConj H m g : cyclotome ↥H m ≃* cyclotome ↥H m`。左作用であること
   (`cyclotomeConj_one` / `cyclotomeConj_mul`)。
4. **内部自己同型は自明に作用する**(`topAbelianizationCME_conj_self` /
   `cyclotomeConj_coe_self`)。すなわち作用は `Γ ⧸ H = Gal(L_H/K)` を経由する。
   ★これが「`Aut(L/K)`-同変」という言い方を正当化する内容である。
5. **`α` に沿った移送の同変性**(`cyclotomeEquiv_subgroupMapCME_conj`)——
   `α : Γ ≃ₜ* Γ'` が誘導する `cyclotome ↥H m ≃* cyclotome ↥(α H) m` は
   `g` の作用を `α g` の作用に移す。★経路 Λ の要石はここである。

## 逸脱の記録

**`Aut(ℤ/p^n) ≅ (ℤ/p^n)^×` の標準同一視は作っていない。**
規模測定の節点記述はこれを Λ2 の内容に数えていたが、下流
(`Found/PGC/CyclotomeTransport.lean`)は `Λ_n(H)` を `ZMod (p^n)` に同一視せず
`μ_{p^n} ⊆ K̄` と直接比較して円分指標を読む設計にしたので、
`ZMod` 経由の中継が不要になった。必要になったら mathlib の
`ZMod.AddAutEquivUnits : AddAut (ZMod n) ≃+ Additive (ZMod n)ˣ` がそのまま使える
(2026-09-06 に在庫確認済み)。
-/

namespace ABC3.Found.PGC

open scoped commutatorElement

/-! ## 1. 円分子の定義 -/

/-- **`A^{ab}` の `m` 捩れ部分群**。`Λ_n(H)` は `cyclotome ↥H (p^n)` として得る。

★引数は位相群 `A` だけである。「`Γ` の部分群 `H`」はデータに入っていない
(モジュール docstring の設計を参照)。 -/
def cyclotome (A : Type*) [Group A] [TopologicalSpace A] [IsTopologicalGroup A] (m : ℕ) :
    Subgroup (TopologicalAbelianization A) :=
  (powMonoidHom m : TopologicalAbelianization A →* TopologicalAbelianization A).ker

@[simp] theorem mem_cyclotome {A : Type*} [Group A] [TopologicalSpace A] [IsTopologicalGroup A]
    {m : ℕ} {x : TopologicalAbelianization A} : x ∈ cyclotome A m ↔ x ^ m = 1 := Iff.rfl

/-! ## 2. 移送 -/

/-- 群同型は `m` 乗写像の核を核に写す。 -/
theorem map_powKer_of_mulEquiv {A B : Type*} [CommGroup A] [CommGroup B] (e : A ≃* B) (m : ℕ) :
    ((powMonoidHom m : A →* A).ker).map e.toMonoidHom = (powMonoidHom m : B →* B).ker := by
  ext y
  refine ⟨?_, ?_⟩
  · rintro ⟨x, hx, rfl⟩
    have hx' : x ^ m = 1 := hx
    show (e x) ^ m = 1
    rw [← map_pow, hx', map_one]
  · intro hy
    refine ⟨e.symm y, ?_, by simp⟩
    show (e.symm y) ^ m = 1
    rw [← map_pow]
    have hy' : y ^ m = 1 := hy
    rw [hy', map_one]

/-- **位相群の同型が誘導する円分子の同型**。 -/
noncomputable def cyclotomeEquiv {A B : Type*} [Group A] [Group B] [TopologicalSpace A]
    [TopologicalSpace B] [IsTopologicalGroup A] [IsTopologicalGroup B] (α : A ≃ₜ* B) (m : ℕ) :
    ↥(cyclotome A m) ≃* ↥(cyclotome B m) :=
  (((topAbelianizationCME α).toMulEquiv).subgroupMap (cyclotome A m)).trans
    (MulEquiv.subgroupCongr (map_powKer_of_mulEquiv _ m))

@[simp] theorem cyclotomeEquiv_coe {A B : Type*} [Group A] [Group B] [TopologicalSpace A]
    [TopologicalSpace B] [IsTopologicalGroup A] [IsTopologicalGroup B] (α : A ≃ₜ* B) (m : ℕ)
    (x : ↥(cyclotome A m)) :
    ((cyclotomeEquiv α m x : ↥(cyclotome B m)) : TopologicalAbelianization B)
      = topAbelianizationCME α (x : TopologicalAbelianization A) := rfl

/-- 点ごとに一致する 2 つの同型は同じ誘導写像を与える。 -/
theorem topAbelianizationCME_congr {A B : Type*} [Group A] [Group B] [TopologicalSpace A]
    [TopologicalSpace B] [IsTopologicalGroup A] [IsTopologicalGroup B] {α β : A ≃ₜ* B}
    (h : ∀ a, α a = β a) (x : TopologicalAbelianization A) :
    topAbelianizationCME α x = topAbelianizationCME β x := by
  induction x using QuotientGroup.induction_on with
  | H g =>
    show (QuotientGroup.mk (α g) : TopologicalAbelianization B) = QuotientGroup.mk (β g)
    rw [h g]

theorem cyclotomeEquiv_congr {A B : Type*} [Group A] [Group B] [TopologicalSpace A]
    [TopologicalSpace B] [IsTopologicalGroup A] [IsTopologicalGroup B] {α β : A ≃ₜ* B}
    (h : ∀ a, α a = β a) (m : ℕ) (x : ↥(cyclotome A m)) :
    cyclotomeEquiv α m x = cyclotomeEquiv β m x :=
  Subtype.ext (topAbelianizationCME_congr h _)

theorem cyclotomeEquiv_trans {A B C : Type*} [Group A] [Group B] [Group C]
    [TopologicalSpace A] [TopologicalSpace B] [TopologicalSpace C]
    [IsTopologicalGroup A] [IsTopologicalGroup B] [IsTopologicalGroup C]
    (α : A ≃ₜ* B) (β : B ≃ₜ* C) (m : ℕ) (x : ↥(cyclotome A m)) :
    cyclotomeEquiv β m (cyclotomeEquiv α m x) = cyclotomeEquiv (α.trans β) m x :=
  Subtype.ext (topAbelianizationCME_trans α β _)

theorem cyclotomeEquiv_refl {A : Type*} [Group A] [TopologicalSpace A] [IsTopologicalGroup A]
    (m : ℕ) (x : ↥(cyclotome A m)) : cyclotomeEquiv (ContinuousMulEquiv.refl A) m x = x :=
  Subtype.ext (topAbelianizationCME_refl _)

/-! ## 3. 共役——正規部分群への `Γ` の作用 -/

/-- **`Γ` の正規部分群 `H` への共役作用**、位相群の同型として。

`MulAut.conjNormal` に連続性を付けただけ。連続性は
`x ↦ g * x * g⁻¹` が `Γ` 上連続で `↥H` が誘導位相を持つことから出る。 -/
noncomputable def conjSubgroupCME {Γ : Type*} [Group Γ] [TopologicalSpace Γ] [IsTopologicalGroup Γ]
    {H : Subgroup Γ} [H.Normal] (g : Γ) : ↥H ≃ₜ* ↥H where
  toMulEquiv := MulAut.conjNormal g
  continuous_toFun :=
    continuous_induced_rng.2
      (((continuous_const.mul continuous_subtype_val).mul continuous_const).congr
        (fun x => (MulAut.conjNormal_apply g x).symm))
  continuous_invFun :=
    continuous_induced_rng.2
      (((continuous_const.mul continuous_subtype_val).mul continuous_const).congr
        (fun x => (MulAut.conjNormal_symm_apply (H := H) g x).symm))

@[simp] theorem conjSubgroupCME_coe {Γ : Type*} [Group Γ] [TopologicalSpace Γ]
    [IsTopologicalGroup Γ] {H : Subgroup Γ} [H.Normal] (g : Γ) (x : ↥H) :
    ((conjSubgroupCME g x : ↥H) : Γ) = g * (x : Γ) * g⁻¹ := MulAut.conjNormal_apply g x

theorem conjSubgroupCME_one {Γ : Type*} [Group Γ] [TopologicalSpace Γ] [IsTopologicalGroup Γ]
    {H : Subgroup Γ} [H.Normal] (x : ↥H) : conjSubgroupCME (H := H) 1 x = x := by
  apply Subtype.ext
  rw [conjSubgroupCME_coe]
  group

theorem conjSubgroupCME_mul {Γ : Type*} [Group Γ] [TopologicalSpace Γ] [IsTopologicalGroup Γ]
    {H : Subgroup Γ} [H.Normal] (g₁ g₂ : Γ) (x : ↥H) :
    conjSubgroupCME (H := H) (g₁ * g₂) x = conjSubgroupCME g₁ (conjSubgroupCME g₂ x) := by
  apply Subtype.ext
  rw [conjSubgroupCME_coe, conjSubgroupCME_coe, conjSubgroupCME_coe]
  group

/-- **内部自己同型はアーベル化に自明に作用する**。

`h ∈ H` による共役は `H^{ab}` の上で恒等。したがって `Γ` の `H^{ab}` への作用は
`Γ ⧸ H = Gal(L_H/K)` を経由する。★これが「`Aut(L/K)`-同変」という言い方の中身。 -/
theorem topAbelianizationCME_conj_self {Γ : Type*} [Group Γ] [TopologicalSpace Γ]
    [IsTopologicalGroup Γ] {H : Subgroup Γ} [H.Normal] (h : ↥H)
    (x : TopologicalAbelianization ↥H) :
    topAbelianizationCME (conjSubgroupCME (H := H) (h : Γ)) x = x := by
  induction x using QuotientGroup.induction_on with
  | H y =>
    show (QuotientGroup.mk (conjSubgroupCME (H := H) (h : Γ) y) : TopologicalAbelianization ↥H)
      = QuotientGroup.mk y
    have hy : conjSubgroupCME (H := H) (h : Γ) y = h * y * h⁻¹ := by
      apply Subtype.ext
      show (h : Γ) * (y : Γ) * (h : Γ)⁻¹ = ((h * y * h⁻¹ : ↥H) : Γ)
      push_cast
      ring_nf
    rw [hy]
    apply QuotientGroup.eq.mpr
    have hc : (h * y * h⁻¹)⁻¹ * y = ⁅h, y⁻¹⁆ := by group
    rw [hc]
    exact Subgroup.le_topologicalClosure _
      (Subgroup.commutator_mem_commutator (Subgroup.mem_top h) (Subgroup.mem_top y⁻¹))

/-! ## 4. 円分子への共役作用 -/

/-- **`Λ_n(H)` への `Γ` の作用**——`g : Γ` による共役が誘導する円分子の自己同型。 -/
noncomputable def cyclotomeConj {Γ : Type*} [Group Γ] [TopologicalSpace Γ] [IsTopologicalGroup Γ]
    (H : Subgroup Γ) [H.Normal] (m : ℕ) (g : Γ) :
    ↥(cyclotome ↥H m) ≃* ↥(cyclotome ↥H m) :=
  cyclotomeEquiv (conjSubgroupCME (H := H) g) m

@[simp] theorem cyclotomeConj_coe {Γ : Type*} [Group Γ] [TopologicalSpace Γ] [IsTopologicalGroup Γ]
    (H : Subgroup Γ) [H.Normal] (m : ℕ) (g : Γ) (x : ↥(cyclotome ↥H m)) :
    ((cyclotomeConj H m g x : ↥(cyclotome ↥H m)) : TopologicalAbelianization ↥H)
      = topAbelianizationCME (conjSubgroupCME (H := H) g) (x : TopologicalAbelianization ↥H) := rfl

theorem cyclotomeConj_one {Γ : Type*} [Group Γ] [TopologicalSpace Γ] [IsTopologicalGroup Γ]
    (H : Subgroup Γ) [H.Normal] (m : ℕ) (x : ↥(cyclotome ↥H m)) :
    cyclotomeConj H m 1 x = x := by
  refine Subtype.ext ?_
  rw [cyclotomeConj_coe,
    topAbelianizationCME_congr (β := ContinuousMulEquiv.refl ↥H) conjSubgroupCME_one]
  exact topAbelianizationCME_refl _

theorem cyclotomeConj_mul {Γ : Type*} [Group Γ] [TopologicalSpace Γ] [IsTopologicalGroup Γ]
    (H : Subgroup Γ) [H.Normal] (m : ℕ) (g₁ g₂ : Γ) (x : ↥(cyclotome ↥H m)) :
    cyclotomeConj H m (g₁ * g₂) x = cyclotomeConj H m g₁ (cyclotomeConj H m g₂ x) := by
  rw [cyclotomeConj, cyclotomeConj, cyclotomeConj, cyclotomeEquiv_trans]
  exact cyclotomeEquiv_congr (fun a => conjSubgroupCME_mul g₁ g₂ a) m x

/-- `H` の元による作用は自明——作用は `Γ ⧸ H` を経由する。 -/
theorem cyclotomeConj_coe_self {Γ : Type*} [Group Γ] [TopologicalSpace Γ] [IsTopologicalGroup Γ]
    (H : Subgroup Γ) [H.Normal] (m : ℕ) (h : ↥H) (x : ↥(cyclotome ↥H m)) :
    cyclotomeConj H m (h : Γ) x = x :=
  Subtype.ext (topAbelianizationCME_conj_self h _)

/-! ## 5. ★同変性——経路 Λ の要石 -/

/-- **`α` による共役の移送**——`α : Γ ≃ₜ* Γ'` の `H` への制限は
`g` による共役を `α g` による共役に移す。 -/
theorem subgroupMapCME_conjSubgroupCME {Γ Γ' : Type*} [Group Γ] [Group Γ']
    [TopologicalSpace Γ] [TopologicalSpace Γ'] [IsTopologicalGroup Γ] [IsTopologicalGroup Γ']
    (α : Γ ≃ₜ* Γ') (H : Subgroup Γ) [H.Normal]
    (_hN : (H.map α.toMulEquiv.toMonoidHom).Normal) (g : Γ) (x : ↥H) :
    subgroupMapCME α H (conjSubgroupCME g x)
      = conjSubgroupCME (H := H.map α.toMulEquiv.toMonoidHom) (α g) (subgroupMapCME α H x) := by
  apply Subtype.ext
  show α ((conjSubgroupCME g x : ↥H) : Γ) = _
  rw [conjSubgroupCME_coe]
  show α (g * (x : Γ) * g⁻¹)
    = α g * (((subgroupMapCME α H x : ↥(H.map α.toMulEquiv.toMonoidHom)) : Γ')) * (α g)⁻¹
  show α (g * (x : Γ) * g⁻¹) = α g * (α (x : Γ)) * (α g)⁻¹
  rw [map_mul, map_mul, map_inv]

/-- **★★★★★★★★★★★★円分子の移送は `Γ` の作用と同変**。

`α : Γ ≃ₜ* Γ'`、`H ⊴ Γ`、`H' := α(H)` について
`cyclotomeEquiv (subgroupMapCME α H) : Λ_m(H) ≃* Λ_m(H')` は
`g` の作用を `α g` の作用に移す:

  `e (g · x) = (α g) · (e x)`。

★`H` はここでは仮定の側にしか現れない。下流(`CyclotomeTransport.lean`)は
これを `∃ H` の内側で使うので、結論の statement に `H` は出てこない。 -/
theorem cyclotomeEquiv_subgroupMapCME_conj {Γ Γ' : Type*} [Group Γ] [Group Γ']
    [TopologicalSpace Γ] [TopologicalSpace Γ'] [IsTopologicalGroup Γ] [IsTopologicalGroup Γ']
    (α : Γ ≃ₜ* Γ') (H : Subgroup Γ) [H.Normal]
    (hN : (H.map α.toMulEquiv.toMonoidHom).Normal) (m : ℕ) (g : Γ)
    (x : ↥(cyclotome ↥H m)) :
    cyclotomeEquiv (subgroupMapCME α H) m (cyclotomeConj H m g x)
      = @cyclotomeConj Γ' _ _ _ (H.map α.toMulEquiv.toMonoidHom) hN m (α g)
          (cyclotomeEquiv (subgroupMapCME α H) m x) := by
  rw [cyclotomeConj, cyclotomeConj, cyclotomeEquiv_trans, cyclotomeEquiv_trans]
  exact cyclotomeEquiv_congr (fun a => subgroupMapCME_conjSubgroupCME α H hN g a) m x

end ABC3.Found.PGC
