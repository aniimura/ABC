import ABC3.Found.PGC.RamificationFiltrationBuild
import ABC3.Found.PGC.AbelianSplitUnramified

/-!
# 不分岐底変換に沿った上付き番号付けの輸送 —— ★★`StageFiltration.compat` が閉じた

`Found/PGC/RamificationFiltrationBuild.lean`(Y19d)は
`Interface/PGC/LocalFieldData.lean` の `RamificationFiltration p` を
**`compat` ただ 1 つ**に還元した(`ramificationFiltrationOfCompat`)。
★★**本ファイルはその `compat` を無条件で埋める。**

到達点は

```
ramificationFiltration (p : ℕ) [Fact p.Prime] : RamificationFiltration p
```

であり、**仮定は 1 つも無い**。★`Skeleton/PGC/Section2.lean` の
`prop_2_1` / `prop_2_2` はこれをそのまま入力に取れる。

## 何が入ったか

`L := K(x) ⊆ L′ := K(x′)`、`G := I(L/K)`、`G′ := I(L′/K)`(絶対惰性群 `I_K` の像)、
`f : G′ ↠ G` を制限射、`H := ker f`、`C := (𝒪_{L′})^H` と置く。

* §1 **抽象核(純群論)** —— 分岐・付値・Galois の語彙が 1 つも出てこない:
  - `phiOf_comp_of_surjective` : ★**Herbrand 関数 `φ` は全射準同型に沿って不変**。
    中身は「ファイバーの大きさが一定」+ `|G₁| = |ker| · |G₂|` だけである。
  - `descendHom` / `descendHom_apply` : `ker θ₁ ≤ ker θ₂` なら `θ₂` は `θ₁` を経由する。
  - `coe_comap_mul_coe_ker` : ★★**`compat` の骨格**。`ψ = δ ∘ ψ′` のとき
    `(ψ′⁻¹ S′) · ker ψ = ψ⁻¹ (δ S′)`。
* §2 **抽象核(輸送)** —— `ramIndex` が対応すれば上付き分岐群も対応する:
  `upperRamificationGroup_eq_comap` / `map_upperRamificationGroup_eq`。
  ★本体の見当「輸送は同型に沿って写すだけ」は**半分当たった**: 実際に要るのは
  同型ではなく**全射準同型**であり(`G′ ↠ G` は同型ではない)、
  `ramIndex` さえ対応すれば `φ`・`ψ`・`G^m` がすべて対応する。
* §3 **抽象核** —— 同変かつ「不分岐」(素元が素元に行く)な環準同型は `ramIndex` を保つ:
  `addVal_ringHom_eq_of_irreducible` / `ramIndex_ringHom_eq`。
* §4–§6 **具体層** —— `restrictGalHom` / `inertiaRestrict` / `stageBaseHom`。
* §7–§8 **`compat`** —— `absGalStage_compat` と
  `ramificationFiltrationOfUnramifiedBaseChange`(不分岐性を仮定に取る形)。
* §9 **`addVal` と `Ideal.ramificationIdx` の橋** ——
  `addVal_algebraMap_eq_ramificationIndex`。★`QuotientRamIndexAverage.lean` が
  「mathlib の `Ideal.ramificationIdx` は `addVal` と繋がっていない」と書いた穴を塞いだ。
  中身は `sSup {n | n ≤ v} = v`(`csSup_Iic`)だけである。
* §10 ★★★**`e = |I|`**(`ramificationIndex_eq_card_inertiaGal`)——
  分岐指数は惰性群の位数に等しい。決め手は
  **`finrank K L₀ = f(L/K)`**(`finrank_fixedField_inertiaGal`):
  `L₀ = K^ur ⊓ L`(Y19b+c)を原始元で `K(w)` と書き、
  `AbelianSplitUnramified.lean` の `exists_unramified_subextension`
  (最大不分岐部分拡大の次数が `f`)と `inertiaDegree_le_of_adjoin_le` で挟む。
* §11 ★★★★★**`unramifiedBaseChange`**(仮定が定理になった)と
  **`ramificationFiltration p`**(仮定なしの `RamificationFiltration p`)。

## ★何が入っていないか(正直な線引き)

★★**`compat` は無条件で閉じた。** `sorry` は 1 つも無く、
`#print axioms ramificationFiltration` は `[propext, Classical.choice, Quot.sound]` である。

★★**近似(`trivialStageFiltration` / `inertiaStageFiltration`)では埋めていない。**
本ファイルが作るのは Y19d の `absGalStage`(本物)であり、近似 2 つと値が食い違うことを
`exists_stage_ne_inertiaStageFiltration` / `stage_ne_trivialStageFiltration` で確認してある。
さらに `coe_ramificationFiltration_mul_coe` により
**`Γ_K^v` は各有限段 `Gal(L/L₀)^v` へ全射する**(退化した族ではこの式は成り立たない)。

★**証明していないこと**(下流に影響しない範囲):
`G^m` が一意化元 `stageUniformizer` と生成元 `StageGenerator` の選択に依らないこと
(Y19d の逸脱 1・2 をそのまま引き継ぐ)。`compat` の証明では両側とも
`(nonempty_stageGenerator K _).some` という**同じ選択**を使うので要らなかった。

## 逸脱の記録(CLAUDE.md「逸脱」)

1. ★`UnramifiedBaseChangeHyp` は「すべての `K` と `K(x) ⊆ K(x′)` の対」について
   述べている(`compat` に要るのは `StageGenerator` が選んだ生成元の対だけ)。
   ★**定理として証明したので仮定の強さは問題にならない。**
2. ★§10 の `isGalois_adjoin` は「`K.closure/K` が Galois」から分離性を降ろしている
   (標数 0 を直接は使っていない)。
3. ★段の一意化元 `stageUniformizer` と生成元 `StageGenerator` の選択は Y19d のものを
   そのまま使う。**`G^m` がその選択に依らないことは本ファイルでも証明していない**
   (Y19d 逸脱 1・2 を引き継ぐ)。★`compat` の証明では両側とも
   `(nonempty_stageGenerator K _).some` という**同じ選択**を使うので、
   生成元の独立性は要らなかった。
4. ★本ファイルは原典が名前を付けていない組み立てノードなので `.src` を持たない
   (`InertiaReduction.lean` / `RamificationFiltrationBuild.lean` と同じ理由)。
   使っている原典項目は Yoshida Corollary 6.13(i)(= Y18 `upperRamification_coe_mul_coe_eq`)
   と Definition 6.12 であり、`.src` はそれぞれの定義元にある。
5. ★**既存の `Found/PGC/*.lean` は 1 行も書き換えていない**(import のみ)。

## 配管の記録(`lean-idioms.md` 行き)

* ★**中間体の 2 層の塔を 1 つも作っていない**(#59 回避)。`Gal(L′/K) → Gal(L/K)` は
  `↥L′` の上に `↥L` を載せるのではなく、`Γ_K` からの 2 本の全射を
  `descendHom` で繋いで作った(`restrictGalHom`)。
* ★**`adjoinIntegersIncl` の係数の `rfl` は 2 層にすると kernel が落ちる**(実測)。
  `val_adjoinIntegersIncl`(1 層・明示の無名構成子)を先に `rfl` で作り、
  `congrArg Subtype.val` で 2 層目に上げると通る。
* ★★**`ABC3.Found.PGC.ker_restrictNormalHom_eq_fixingSubgroup` は同名が 2 つある**
  (`AbsGalRamificationFiltration.lean:492` と `LubinTateClosure.lean:127`)。
  両方が import 圏に入ると曖昧になって落ちる。★mathlib の
  `IntermediateField.restrictNormalHom_ker` を使えば回避できる。
* ★`ℕ∞` の切り詰め引き算・除算は 1 箇所も書いていない(#102)。
* ★`set M := …` のあとに `α : ↥M` を取ると `rw [hM]` の motive が壊れる(#135)。
  `∃ w : K.closure, adjoin K {w} = M` の形で **`↥M` に依らない `w` を先に取り出す**と通る。
-/

namespace ABC3.Found.PGC

open Finset IsLocalRing IsDiscreteValuationRing
open scoped Pointwise

/-! ## §1 抽象核 —— 純群論(分岐・付値・Galois の語彙は 1 つも出てこない) -/

section Pure

variable {G₁ G₂ : Type*} [Group G₁] [Group G₂] [Fintype G₁] [Fintype G₂]

omit [Fintype G₂] in
/-- 全射準同型のファイバーはどれも `|ker|` 個。 -/
theorem card_fiber_monoidHom {f : G₁ →* G₂} (hf : Function.Surjective f) [DecidableEq G₂]
    (τ : G₂) : #{σ : G₁ | f σ = τ} = Nat.card f.ker := by
  classical
  rw [MonoidHom.card_fiber_eq_of_mem_range f (hf τ) ⟨1, map_one f⟩,
    Nat.card_eq_fintype_card, Fintype.card_subtype]
  simp [MonoidHom.mem_ker]

/-- 引き戻した関数の和は `|ker|` 倍。 -/
theorem sum_comp_monoidHom_surjective {M : Type*} [AddCommMonoid M]
    {f : G₁ →* G₂} (hf : Function.Surjective f) (F : G₂ → M) :
    ∑ σ : G₁, F (f σ) = ∑ τ : G₂, (Nat.card f.ker) • F τ := by
  classical
  rw [← Finset.sum_fiberwise_of_maps_to (g := f) (t := (univ : Finset G₂))
    (fun x _ => mem_univ (f x))]
  refine Finset.sum_congr rfl (fun τ _ => ?_)
  rw [Finset.sum_congr rfl (fun σ hσ => congrArg F (mem_filter.1 hσ).2), Finset.sum_const,
    card_fiber_monoidHom hf τ]

omit [Fintype G₁] [Fintype G₂] in
/-- `|G₁| = |ker f| · |G₂|`。 -/
theorem card_eq_card_ker_mul {f : G₁ →* G₂} (hf : Function.Surjective f) :
    Nat.card G₁ = Nat.card f.ker * Nat.card G₂ := by
  classical
  have h := Subgroup.card_mul_index f.ker
  have h2 : f.ker.index = Nat.card G₂ := by
    rw [← Nat.card_congr (QuotientGroup.quotientKerEquivOfSurjective f hf).toEquiv]
    rfl
  rw [← h, h2]

/-- ★★★**抽象核** —— Herbrand 関数の器 `φ_f` は**全射準同型に沿って不変**である。

`φ_f(n) = −1 + (1/|S|) Σ_τ min{f(τ), n+1}` の分子が `|ker|` 倍、分母も `|ker|` 倍に
なるだけである。★これが「上付き番号付けは商と両立する」ことの数え上げの中身であり、
分岐・付値・Galois の語彙は 1 つも出てこない。 -/
theorem phiOf_comp_of_surjective {f : G₁ →* G₂} (hf : Function.Surjective f) (g : G₂ → ℕ∞)
    (n : ℝ) : phiOf (fun σ : G₁ => g (f σ)) n = phiOf g n := by
  have hk : (0 : ℕ) < Nat.card f.ker := Nat.card_pos
  have hk' : (Nat.card f.ker : ℝ) ≠ 0 := Nat.cast_ne_zero.2 hk.ne'
  have hsum : ∑ σ : G₁, truncENat (g (f σ)) (n + 1)
      = (Nat.card f.ker : ℝ) * ∑ τ : G₂, truncENat (g τ) (n + 1) := by
    rw [sum_comp_monoidHom_surjective hf (fun τ => truncENat (g τ) (n + 1)), ← Finset.smul_sum,
      nsmul_eq_mul]
  rw [phiOf, phiOf, hsum, card_eq_card_ker_mul hf, Nat.cast_mul, mul_div_mul_left _ _ hk']

end Pure

section Descend

variable {Γ G₁ G₂ : Type*} [Group Γ] [Group G₁] [Group G₂]

/-- ★★**抽象核** —— `θ₁` が全射で `ker θ₁ ≤ ker θ₂` なら `θ₂` は `G₁` を経由する。

★★これで**中間体の 2 層の塔を作らずに** `Gal(L′/K) → Gal(L/K)` が作れる
(`lean-idioms.md` #59 回避)。 -/
noncomputable def descendHom (θ₁ : Γ →* G₁) (h₁ : Function.Surjective θ₁) (θ₂ : Γ →* G₂)
    (hker : θ₁.ker ≤ θ₂.ker) : G₁ →* G₂ :=
  (QuotientGroup.lift θ₁.ker θ₂ (fun _ hg => hker hg)).comp
    (QuotientGroup.quotientKerEquivOfSurjective θ₁ h₁).symm.toMonoidHom

@[simp] theorem descendHom_apply (θ₁ : Γ →* G₁) (h₁ : Function.Surjective θ₁) (θ₂ : Γ →* G₂)
    (hker : θ₁.ker ≤ θ₂.ker) (g : Γ) : descendHom θ₁ h₁ θ₂ hker (θ₁ g) = θ₂ g := by
  have h : (QuotientGroup.quotientKerEquivOfSurjective θ₁ h₁).symm (θ₁ g)
      = (QuotientGroup.mk g : Γ ⧸ θ₁.ker) :=
    (MulEquiv.symm_apply_eq _).2 rfl
  rw [descendHom, MonoidHom.comp_apply, MulEquiv.coe_toMonoidHom, h]
  rfl

theorem surjective_descendHom (θ₁ : Γ →* G₁) (h₁ : Function.Surjective θ₁) (θ₂ : Γ →* G₂)
    (hker : θ₁.ker ≤ θ₂.ker) (h₂ : Function.Surjective θ₂) :
    Function.Surjective (descendHom θ₁ h₁ θ₂ hker) := by
  intro y
  obtain ⟨g, rfl⟩ := h₂ y
  exact ⟨θ₁ g, descendHom_apply θ₁ h₁ θ₂ hker g⟩

end Descend

section CompatCore

variable {Γ Q' Q : Type*} [Group Γ] [Group Q'] [Group Q]

/-- ★★★**抽象核 —— `StageFiltration.compat` の骨格**。

2 つの全射 `ψ′ : Γ ↠ Q′`、`ψ : Γ ↠ Q` が `ψ = δ ∘ ψ′` で繋がっているとき、
`Q′` の部分群 `S′` の引き戻しに `ker ψ` を掛けると `δ(S′)` の引き戻しになる:

`(ψ′⁻¹ S′) · ker ψ = ψ⁻¹ (δ S′)`。

★分岐・付値・Galois の語彙は 1 つも出てこない。★`Γ` の有限性も位相も要らない。 -/
theorem coe_comap_mul_coe_ker {ψ' : Γ →* Q'} {ψ : Γ →* Q} {δ : Q' →* Q}
    (hψ' : Function.Surjective ψ') (hcomp : ∀ γ : Γ, ψ γ = δ (ψ' γ)) (S' : Subgroup Q') :
    ((Subgroup.comap ψ' S' : Subgroup Γ) : Set Γ) * ((ψ.ker : Subgroup Γ) : Set Γ)
      = ((Subgroup.comap ψ (Subgroup.map δ S') : Subgroup Γ) : Set Γ) := by
  rw [← Subgroup.mul_normal]
  congr 1
  refine le_antisymm (sup_le (fun γ hγ => ?_) (fun γ hγ => ?_)) (fun γ hγ => ?_)
  · exact Subgroup.mem_comap.2 (by rw [hcomp γ]; exact ⟨ψ' γ, Subgroup.mem_comap.1 hγ, rfl⟩)
  · exact Subgroup.mem_comap.2 (by rw [MonoidHom.mem_ker.1 hγ]; exact one_mem _)
  · obtain ⟨s, hs, hδ⟩ := Subgroup.mem_comap.1 hγ
    obtain ⟨γ₀, rfl⟩ := hψ' s
    have hker : γ * γ₀⁻¹ ∈ ψ.ker := by
      have heq : ψ γ₀ = ψ γ := by rw [hcomp γ₀, hδ]
      rw [MonoidHom.mem_ker, map_mul, map_inv, heq, mul_inv_cancel]
    have h₀ : γ₀ ∈ Subgroup.comap ψ' S' := Subgroup.mem_comap.2 hs
    have hsplit : γ = (γ * γ₀⁻¹) * γ₀ := by group
    rw [hsplit]
    exact mul_mem (Subgroup.mem_sup_right hker) (Subgroup.mem_sup_left h₀)

end CompatCore

/-! ## §2 抽象核 —— `ramIndex` が対応すれば上付き分岐群も対応する

★★**本体の見当「輸送は同型に沿って写すだけ」の訂正**: 実際に要るのは**同型ではなく
全射準同型**である(`G′ ↠ G` は同型ではない)。`ramIndex` さえ対応すれば
`G_n`(実数添字)・`φ`・`ψ`・`G^m` がこの順で全部対応する。 -/

section Transport

variable {B₁ B₂ : Type*} [CommRing B₁] [IsDomain B₁] [IsDiscreteValuationRing B₁]
  [CommRing B₂] [IsDomain B₂] [IsDiscreteValuationRing B₂]
variable {G₁ G₂ : Type*} [Group G₁] [Group G₂] [Fintype G₁] [Fintype G₂]
  [MulSemiringAction G₁ B₁] [MulSemiringAction G₂ B₂]

omit [Fintype G₁] [Fintype G₂] in
/-- 実数添字の下付き分岐群は引き戻しで対応する(全射性も要らない)。 -/
theorem ramificationGroupReal_eq_comap (f : G₁ →* G₂) {α₁ : B₁} {α₂ : B₂}
    (h : ∀ σ : G₁, ramIndex α₁ σ = ramIndex α₂ (f σ)) (n : ℝ) :
    ramificationGroupReal α₁ n = Subgroup.comap f (ramificationGroupReal α₂ n) := by
  ext σ
  simp only [Subgroup.mem_comap, mem_ramificationGroupReal, h σ]

/-- Herbrand 関数 `φ_G` が一致する(§1 の抽象核に `ramIndex` を代入するだけ)。 -/
theorem herbrandPhiGroup_eq_of_ramIndex {f : G₁ →* G₂} (hf : Function.Surjective f)
    {α₁ : B₁} {α₂ : B₂} (h : ∀ σ : G₁, ramIndex α₁ σ = ramIndex α₂ (f σ)) (n : ℝ) :
    herbrandPhiGroup G₁ α₁ n = herbrandPhiGroup G₂ α₂ n := by
  rw [herbrandPhiGroup, herbrandPhiGroup,
    show (fun σ : G₁ => ramIndex α₁ σ) = (fun σ : G₁ => ramIndex α₂ (f σ)) from funext h]
  exact phiOf_comp_of_surjective hf _ n

/-- 逆関数 `ψ_G` も一致する(`Function.invFun` は関数だけで決まるから)。 -/
theorem herbrandPsiGroup_eq_of_ramIndex {f : G₁ →* G₂} (hf : Function.Surjective f)
    {α₁ : B₁} {α₂ : B₂} (h : ∀ σ : G₁, ramIndex α₁ σ = ramIndex α₂ (f σ)) (m : ℝ) :
    herbrandPsiGroup G₁ α₁ m = herbrandPsiGroup G₂ α₂ m := by
  rw [herbrandPsiGroup, herbrandPsiGroup,
    show herbrandPhiGroup G₁ α₁ = herbrandPhiGroup G₂ α₂ from
      funext (herbrandPhiGroup_eq_of_ramIndex hf h)]

/-- ★★★**抽象核の到達点** —— 上付き分岐群は引き戻しで対応する。 -/
theorem upperRamificationGroup_eq_comap {f : G₁ →* G₂} (hf : Function.Surjective f)
    {α₁ : B₁} {α₂ : B₂} (h : ∀ σ : G₁, ramIndex α₁ σ = ramIndex α₂ (f σ)) (m : ℝ) :
    upperRamificationGroup G₁ α₁ m = Subgroup.comap f (upperRamificationGroup G₂ α₂ m) := by
  rw [upperRamificationGroup_def, upperRamificationGroup_def,
    herbrandPsiGroup_eq_of_ramIndex hf h m, ramificationGroupReal_eq_comap f h]

/-- ★★★像の形(`compat` で使うのはこちら)。 -/
theorem map_upperRamificationGroup_eq {f : G₁ →* G₂} (hf : Function.Surjective f)
    {α₁ : B₁} {α₂ : B₂} (h : ∀ σ : G₁, ramIndex α₁ σ = ramIndex α₂ (f σ)) (m : ℝ) :
    Subgroup.map f (upperRamificationGroup G₁ α₁ m) = upperRamificationGroup G₂ α₂ m := by
  rw [upperRamificationGroup_eq_comap hf h m, Subgroup.map_comap_eq_self_of_surjective hf]

end Transport

/-! ## §3 抽象核 —— 同変な「不分岐」環準同型は `ramIndex` を保つ -/

section RamIndexMap

variable {A B : Type*} [CommRing A] [IsDomain A] [IsDiscreteValuationRing A]
  [CommRing B] [IsDomain B] [IsDiscreteValuationRing B]

/-- ★★**「不分岐」の使いやすい形** —— `A` の素元が `B` でも素元なら `addVal` は保たれる。

★★これが「不分岐底変換」の本体である(`e = 1` を素元の言葉で書いただけ)。 -/
theorem addVal_ringHom_eq_of_irreducible (φ : A →+* B) {ϖ : A} (hϖ : Irreducible ϖ)
    (h : Irreducible (φ ϖ)) (a : A) : addVal B (φ a) = addVal A a := by
  letI := φ.toAlgebra
  have h1 : addVal B (algebraMap A B ϖ) = 1 := addVal_uniformizer h
  show addVal B (algebraMap A B a) = addVal A a
  rw [addVal_map_eq_mul_addVal hϖ (by rw [h1]; exact one_ne_zero) a, h1, one_mul]

variable {G₁ G₂ : Type*} [Group G₁] [Group G₂] [MulSemiringAction G₁ B] [MulSemiringAction G₂ A]

/-- ★★**`ramIndex` の輸送** —— 同変かつ `addVal` を保つ環準同型に沿って `i(σ)` は等しい。

★`ramIndex α σ = addVal (σ•α − α)` なので、証明は
「`σ • φ a = φ (f σ • a)`」と「`addVal ∘ φ = addVal`」を順に当てるだけ(3 語)。 -/
theorem ramIndex_ringHom_eq (φ : A →+* B) (f : G₁ →* G₂)
    (hval : ∀ a : A, addVal B (φ a) = addVal A a)
    (hequiv : ∀ (σ : G₁) (a : A), σ • φ a = φ (f σ • a))
    (α : A) (σ : G₁) : ramIndex (φ α) σ = ramIndex α (f σ) := by
  rw [ramIndex, ramIndex, hequiv, ← map_sub, hval]

end RamIndexMap

/-! ## §4 具体層(1) —— 中間体の対に沿った制限射

★★**塔を作らない**。`Gal(L′/K) → Gal(L/K)` は `Γ_K` からの 2 本の全射を
§1 の `descendHom` で繋いで作る(`lean-idioms.md` #59 回避)。 -/

open ABC3.Skeleton.PGC ABC3.Interface.PGC
open scoped NNReal Valued

attribute [local instance] isDiscreteValuationRing_adjoinIntegers

variable {p : ℕ} [Fact p.Prime]

/-- 制限射の値は大きい体の上では元の自己同型の値である。 -/
theorem coe_restrictNormalHom_apply {F E : Type*} [Field F] [Field E] [Algebra F E]
    (L : IntermediateField F E) [Normal F L] (σ : E ≃ₐ[F] E) (w : ↥L) :
    ((AlgEquiv.restrictNormalHom (F := F) (K₁ := E) (L : Type _) σ w : ↥L) : E) = σ (w : E) := by
  have hc := AlgEquiv.restrictNormal_commutes σ (L : Type _) w
  simp only [AlgEquiv.restrictNormalHom, MonoidHom.mk'_apply]
  exact hc

section RestrictGal

variable (K : PAdicLocalField p) {L L' : IntermediateField K.carrier K.closure}
  [Normal K.carrier L] [Normal K.carrier L']

/-- ★★`L ⊆ L′` に沿った制限射 `Gal(L′/K) → Gal(L/K)`。 -/
noncomputable def restrictGalHom (hle : L ≤ L') :
    (L' ≃ₐ[K.carrier] L') →* (L ≃ₐ[K.carrier] L) :=
  descendHom (AlgEquiv.restrictNormalHom (F := K.carrier) (K₁ := K.closure) (L' : Type _))
    (by haveI := isGalois_closure K; exact surjective_restrictNormalHom L')
    (AlgEquiv.restrictNormalHom (F := K.carrier) (K₁ := K.closure) (L : Type _))
    (by rw [IntermediateField.restrictNormalHom_ker (K := K.carrier) (L := K.closure) (E := L'),
          IntermediateField.restrictNormalHom_ker (K := K.carrier) (L := K.closure) (E := L)]
        exact IntermediateField.fixingSubgroup_le hle)

@[simp] theorem restrictGalHom_apply (hle : L ≤ L') (γ : K.absGal) :
    restrictGalHom K hle
        (AlgEquiv.restrictNormalHom (F := K.carrier) (K₁ := K.closure) (L' : Type _) γ)
      = AlgEquiv.restrictNormalHom (F := K.carrier) (K₁ := K.closure) (L : Type _) γ :=
  descendHom_apply _ _ _ _ γ

theorem comp_restrictGalHom (hle : L ≤ L') :
    (restrictGalHom K hle).comp
        (AlgEquiv.restrictNormalHom (F := K.carrier) (K₁ := K.closure) (L' : Type _))
      = AlgEquiv.restrictNormalHom (F := K.carrier) (K₁ := K.closure) (L : Type _) :=
  MonoidHom.ext (restrictGalHom_apply K hle)

/-- ★惰性群は惰性群の上へ写る —— どちらも同じ `I_K` の像だから。 -/
theorem map_restrictGalHom_inertiaGal (hle : L ≤ L') :
    Subgroup.map (restrictGalHom K hle) (inertiaGal K L') = inertiaGal K L := by
  rw [inertiaGal, inertiaGal, Subgroup.map_map, comp_restrictGalHom]

/-- ★★惰性群の間の制限射 `I(L′/K) ↠ I(L/K)`。★★**その核が Corollary 6.13(i) の `H`** である。 -/
noncomputable def inertiaRestrict (hle : L ≤ L') :
    ↥(inertiaGal K L') →* ↥(inertiaGal K L) :=
  MonoidHom.codRestrict ((restrictGalHom K hle).comp (inertiaGal K L').subtype) _
    (fun σ => by
      rw [← map_restrictGalHom_inertiaGal K hle]
      exact ⟨(σ : L' ≃ₐ[K.carrier] L'), σ.2, rfl⟩)

@[simp] theorem coe_inertiaRestrict (hle : L ≤ L') (σ : ↥(inertiaGal K L')) :
    ((inertiaRestrict K hle σ : ↥(inertiaGal K L)) : L ≃ₐ[K.carrier] L)
      = restrictGalHom K hle (σ : L' ≃ₐ[K.carrier] L') := rfl

theorem surjective_inertiaRestrict (hle : L ≤ L') :
    Function.Surjective (inertiaRestrict K hle) := by
  rintro ⟨τ, hτ⟩
  rw [← map_restrictGalHom_inertiaGal K hle] at hτ
  obtain ⟨σ, hσ, hστ⟩ := hτ
  exact ⟨⟨σ, hσ⟩, Subtype.ext hστ⟩

theorem subtype_comp_inertiaRestrict (hle : L ≤ L') :
    (inertiaGal K L).subtype.comp (inertiaRestrict K hle)
      = (restrictGalHom K hle).comp (inertiaGal K L').subtype := rfl

end RestrictGal

/-! ## §5 具体層(2) —— 整数環の包含は同変

★`adjoinIntegersIncl` の係数を 2 層いっぺんに `rfl` で潰すと **kernel が落ちる**(実測)。
1 層ずつ上げる。 -/

section Incl

variable (K : PAdicLocalField p) {x x' : K.closure}

/-- ★1 層だけの `rfl`(★2 層にすると kernel deterministic timeout になる)。 -/
theorem val_adjoinIntegersIncl (hle : IntermediateField.adjoin K.carrier ({x} : Set K.closure)
      ≤ IntermediateField.adjoin K.carrier ({x'} : Set K.closure)) (z : adjoinIntegers K x) :
    ((adjoinIntegersIncl K hle z).1 :
        IntermediateField.adjoin K.carrier ({x'} : Set K.closure))
      = ⟨(z.1 : K.closure), hle z.1.2⟩ := rfl

/-- 2 層目は `congrArg` で上げる(★`rfl` で書かないこと)。 -/
@[simp] theorem coe_adjoinIntegersIncl
    (hle : IntermediateField.adjoin K.carrier ({x} : Set K.closure)
      ≤ IntermediateField.adjoin K.carrier ({x'} : Set K.closure)) (z : adjoinIntegers K x) :
    (((adjoinIntegersIncl K hle z : adjoinIntegers K x') :
        IntermediateField.adjoin K.carrier ({x'} : Set K.closure)) : K.closure)
      = ((z : IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure) :=
  congrArg Subtype.val (val_adjoinIntegersIncl K hle z)

variable
  [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
  [Normal K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
  [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x'} : Set K.closure))]
  [Normal K.carrier (IntermediateField.adjoin K.carrier ({x'} : Set K.closure))]

/-- `Γ_K` の元で書いた同変性(両辺とも `K.closure` の上では `γ` の値)。 -/
theorem smul_adjoinIntegersIncl_absGal
    (hle : IntermediateField.adjoin K.carrier ({x} : Set K.closure)
      ≤ IntermediateField.adjoin K.carrier ({x'} : Set K.closure)) (γ : K.absGal)
    (z : adjoinIntegers K x) :
    (AlgEquiv.restrictNormalHom (F := K.carrier) (K₁ := K.closure)
        ((IntermediateField.adjoin K.carrier ({x'} : Set K.closure)) : Type _) γ)
        • adjoinIntegersIncl K hle z
      = adjoinIntegersIncl K hle
        ((AlgEquiv.restrictNormalHom (F := K.carrier) (K₁ := K.closure)
          ((IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : Type _) γ) • z) := by
  refine Subtype.ext (Subtype.ext ?_)
  simp only [coe_smul_adjoinIntegers, coe_restrictNormalHom_apply, coe_adjoinIntegersIncl]

/-- ★★**整数環の包含は `restrictGalHom` について同変**。 -/
theorem smul_adjoinIntegersIncl
    (hle : IntermediateField.adjoin K.carrier ({x} : Set K.closure)
      ≤ IntermediateField.adjoin K.carrier ({x'} : Set K.closure))
    (τ : (IntermediateField.adjoin K.carrier ({x'} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x'} : Set K.closure)))
    (z : adjoinIntegers K x) :
    τ • adjoinIntegersIncl K hle z
      = adjoinIntegersIncl K hle (restrictGalHom K hle τ • z) := by
  haveI := isGalois_closure K
  obtain ⟨γ, rfl⟩ := surjective_restrictNormalHom
    (IntermediateField.adjoin K.carrier ({x'} : Set K.closure)) τ
  rw [restrictGalHom_apply]
  exact smul_adjoinIntegersIncl_absGal K hle γ z

end Incl

/-! ## §6 具体層(3) —— `𝒪_L → C` と、残る 1 本の不分岐性 -/

section StageBase

variable (K : PAdicLocalField p) {x x' : K.closure}
  [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
  [Normal K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
  [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x'} : Set K.closure))]
  [Normal K.carrier (IntermediateField.adjoin K.carrier ({x'} : Set K.closure))]

/-- ★★**`𝒪_L → C := (𝒪_{L′})^{ker f}`** —— 像が固定環に入るのは §5 の同変性から即座に出る
(`ker f` の元は `𝒪_L` の上で自明に働く)。★古典的には `C = 𝒪_{L·L′₀}` である。 -/
noncomputable def stageBaseHom (hle : IntermediateField.adjoin K.carrier ({x} : Set K.closure)
      ≤ IntermediateField.adjoin K.carrier ({x'} : Set K.closure)) :
    adjoinIntegers K x →+*
      ↥(fixedRing (adjoinIntegers K x') (inertiaRestrict K hle).ker) :=
  RingHom.codRestrict (adjoinIntegersRingHom K hle) _ (fun z => mem_fixedRing.2 (fun ρ hρ => by
    show ρ.1 • adjoinIntegersIncl K hle z = adjoinIntegersIncl K hle z
    rw [smul_adjoinIntegersIncl K hle ρ.1 z]
    have h1 : restrictGalHom K hle ρ.1 = 1 := by
      rw [← coe_inertiaRestrict K hle ρ, MonoidHom.mem_ker.1 hρ]
      rfl
    rw [h1, one_smul]))

@[simp] theorem coe_stageBaseHom (hle : IntermediateField.adjoin K.carrier ({x} : Set K.closure)
      ≤ IntermediateField.adjoin K.carrier ({x'} : Set K.closure)) (z : adjoinIntegers K x) :
    ((stageBaseHom K hle z : ↥(fixedRing (adjoinIntegers K x') (inertiaRestrict K hle).ker)) :
        adjoinIntegers K x') = adjoinIntegersIncl K hle z := rfl

/-- ★★`stageBaseHom` は `inertiaRestrict` について同変(§3 の `hequiv` を供給する)。 -/
theorem smul_stageBaseHom (hle : IntermediateField.adjoin K.carrier ({x} : Set K.closure)
      ≤ IntermediateField.adjoin K.carrier ({x'} : Set K.closure))
    (σ : ↥(inertiaGalAdjoin K x')) (z : adjoinIntegers K x) :
    σ • stageBaseHom K hle z = stageBaseHom K hle (inertiaRestrict K hle σ • z) := by
  refine Subtype.ext ?_
  show σ.1 • adjoinIntegersIncl K hle z
    = adjoinIntegersIncl K hle ((inertiaRestrict K hle σ).1 • z)
  rw [smul_adjoinIntegersIncl K hle σ.1 z, coe_inertiaRestrict]

end StageBase

/-- ★★★**不分岐底変換の仮定** —— `compat` に残る唯一の数学。

`K(x) ⊆ K(x′)` のとき、`𝒪_{K(x)}` の段の素元 `π` は
`C = (𝒪_{K(x′)})^{ker(I(K(x′)/K) ↠ I(K(x)/K))}` の中でも素元である。

★古典的には `C = 𝒪_{L·L′₀}`(`L′₀ = L′ ∩ K^ur`)であり、`L·L′₀/L` は
**不分岐底変換**なので `e = 1`、すなわち `π` は `C` の素元のままである。
★★これは**近似ではなく仮定**である。証明の分解はファイル冒頭「次のノードのための手順書」。 -/
def UnramifiedBaseChangeHyp (p : ℕ) [Fact p.Prime] : Prop :=
  ∀ (K : PAdicLocalField p) (x x' : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    [Normal K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x'} : Set K.closure))]
    [Normal K.carrier (IntermediateField.adjoin K.carrier ({x'} : Set K.closure))]
    (hle : IntermediateField.adjoin K.carrier ({x} : Set K.closure)
      ≤ IntermediateField.adjoin K.carrier ({x'} : Set K.closure)),
    Irreducible (stageBaseHom K hle (stageUniformizer K x))

/-! ## §7 上付き分岐群の輸送(具体層) -/

section StageTransport

variable (K : PAdicLocalField p) {x x' : K.closure}
  [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
  [Normal K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
  [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x'} : Set K.closure))]
  [Normal K.carrier (IntermediateField.adjoin K.carrier ({x'} : Set K.closure))]

theorem irreducible_stageUniformizer : Irreducible (stageUniformizer K x) :=
  (IsDiscreteValuationRing.irreducible_iff_uniformizer _).2
    (maximalIdeal_eq_span_stageUniformizer K x)

/-- ★★★★**本ファイルの心臓** —— 段の上付き分岐群は制限射で写り合う:

`f (Gal(L′/L′₀)^v) = Gal(L/L₀)^v`。

段取りは 3 つ:
* Y18/Y19d の Corollary 6.13(i)(`stage_upperRamification_coe_mul_coe_eq`)で
  `G′^{π′,v} · H = G′^{ϖ,v}`(`ϖ := stageBaseHom π`、`H := ker f`)。
* `map f` は `H` を潰すので左辺の像は `f (G′^{π′,v})`(`map_sup_of_le_ker`)。
* §2 の抽象核 + §3 の `ramIndex` 輸送で右辺の像が `G^{π,v}`。
  ★ここでだけ不分岐性 `hirr` を使う。 -/
theorem map_inertiaRestrict_upperRamificationGroup
    (hle : IntermediateField.adjoin K.carrier ({x} : Set K.closure)
      ≤ IntermediateField.adjoin K.carrier ({x'} : Set K.closure))
    (hirr : Irreducible (stageBaseHom K hle (stageUniformizer K x))) (v : ℝ) :
    Subgroup.map (inertiaRestrict K hle)
        (upperRamificationGroup ↥(inertiaGalAdjoin K x') (stageUniformizer K x') v)
      = upperRamificationGroup ↥(inertiaGalAdjoin K x) (stageUniformizer K x) v := by
  haveI : Fintype ↥(inertiaRestrict K hle).ker := Fintype.ofFinite _
  have hval : ∀ a : adjoinIntegers K x,
      addVal ↥(fixedRing (adjoinIntegers K x') (inertiaRestrict K hle).ker)
        (stageBaseHom K hle a) = addVal (adjoinIntegers K x) a :=
    fun a => addVal_ringHom_eq_of_irreducible (stageBaseHom K hle)
      (irreducible_stageUniformizer K) hirr a
  have hram : ∀ σ : ↥(inertiaGalAdjoin K x'),
      ramIndex (stageBaseHom K hle (stageUniformizer K x)) σ
        = ramIndex (stageUniformizer K x) (inertiaRestrict K hle σ) :=
    fun σ => ramIndex_ringHom_eq (stageBaseHom K hle) (inertiaRestrict K hle) hval
      (fun τ b => smul_stageBaseHom K hle τ b) _ σ
  have h6 := stage_upperRamification_coe_mul_coe_eq K x' (inertiaRestrict K hle).ker hirr v
  have hsup : upperRamificationGroup ↥(inertiaGalAdjoin K x') (stageUniformizer K x') v
        ⊔ (inertiaRestrict K hle).ker
      = upperRamificationGroup ↥(inertiaGalAdjoin K x')
          (stageBaseHom K hle (stageUniformizer K x)) v :=
    SetLike.coe_injective (by rw [Subgroup.mul_normal]; exact h6)
  rw [← map_sup_of_le_ker (inertiaRestrict K hle) _ _ le_rfl, hsup]
  exact map_upperRamificationGroup_eq (surjective_inertiaRestrict K hle) hram v

/-- `Gal(L/K)` の中で見た形。 -/
theorem map_restrictGalHom_stageUpperRamification
    (hle : IntermediateField.adjoin K.carrier ({x} : Set K.closure)
      ≤ IntermediateField.adjoin K.carrier ({x'} : Set K.closure))
    (hirr : Irreducible (stageBaseHom K hle (stageUniformizer K x))) (v : ℝ) :
    Subgroup.map (restrictGalHom K hle) (stageUpperRamification K x' v)
      = stageUpperRamification K x v := by
  rw [stageUpperRamification, stageUpperRamification, Subgroup.map_map,
    ← subtype_comp_inertiaRestrict K hle, ← Subgroup.map_map,
    map_inertiaRestrict_upperRamificationGroup K hle hirr v]

end StageTransport

/-! ## §8 `compat` と `RamificationFiltration p` -/

/-- `N ≤ M` なら対応する中間体は `L_M ⊆ L_N`(無限次 Galois 対応)。 -/
theorem le_of_stageGenerator_le (K : PAdicLocalField p) {N M : Subgroup K.absGal}
    (gN : StageGenerator K N) (gM : StageGenerator K M) (hNM : N ≤ M) :
    IntermediateField.adjoin K.carrier ({gM.gen} : Set K.closure)
      ≤ IntermediateField.adjoin K.carrier ({gN.gen} : Set K.closure) := by
  haveI := isGalois_closure K
  have h1 : IntermediateField.fixedField M
      = IntermediateField.adjoin K.carrier ({gM.gen} : Set K.closure) := by
    have hf := InfiniteGalois.fixedField_fixingSubgroup
      (IntermediateField.adjoin K.carrier ({gM.gen} : Set K.closure))
    rw [gM.fixingSubgroup_eq] at hf
    exact hf
  have h2 : IntermediateField.fixedField N
      = IntermediateField.adjoin K.carrier ({gN.gen} : Set K.closure) := by
    have hf := InfiniteGalois.fixedField_fixingSubgroup
      (IntermediateField.adjoin K.carrier ({gN.gen} : Set K.closure))
    rw [gN.fixingSubgroup_eq] at hf
    exact hf
  rw [← h1, ← h2]
  intro y hy
  rw [IntermediateField.mem_fixedField_iff] at hy ⊢
  exact fun σ hσ => hy σ (hNM hσ)

/-- ★★★**段データの `compat`**(生成元つきの形)。

★**生成元の取り方の独立性は要らない**: `compat` の両側で
`StageGenerator` を 1 つずつ固定して使うだけである。 -/
theorem stage_mul_coe_eq (h : UnramifiedBaseChangeHyp p) (K : PAdicLocalField p)
    {N M : Subgroup K.absGal} (gN : StageGenerator K N) (gM : StageGenerator K M) (hNM : N ≤ M)
    (v : ℝ) :
    ((gN.stage v : Subgroup K.absGal) : Set K.absGal) * (M : Set K.absGal)
      = ((gM.stage v : Subgroup K.absGal) : Set K.absGal) := by
  letI := gN.finiteDimensional
  letI := gN.normal
  letI := gM.finiteDimensional
  letI := gM.normal
  haveI := isGalois_closure K
  have hle := le_of_stageGenerator_le K gN gM hNM
  have hM : (AlgEquiv.restrictNormalHom (F := K.carrier) (K₁ := K.closure)
      ((IntermediateField.adjoin K.carrier ({gM.gen} : Set K.closure)) : Type _)).ker = M := by
    rw [IntermediateField.restrictNormalHom_ker (K := K.carrier) (L := K.closure),
      gM.fixingSubgroup_eq]
  have key := coe_comap_mul_coe_ker (δ := restrictGalHom K hle)
    (ψ' := AlgEquiv.restrictNormalHom (F := K.carrier) (K₁ := K.closure)
      ((IntermediateField.adjoin K.carrier ({gN.gen} : Set K.closure)) : Type _))
    (by exact surjective_restrictNormalHom _)
    (fun γ => (restrictGalHom_apply K hle γ).symm) (stageUpperRamification K gN.gen v)
  rw [map_restrictGalHom_stageUpperRamification K hle (h K gM.gen gN.gen hle) v, hM] at key
  exact key

/-- ★★★★**`StageFiltration.compat`** —— Y19d が残した唯一の穴。 -/
theorem absGalStage_compat (h : UnramifiedBaseChangeHyp p) (K : PAdicLocalField p)
    ⦃M : Subgroup K.absGal⦄ (hM : M ∈ openNormalBase K.absGal) ⦃N : Subgroup K.absGal⦄
    (hN : N ∈ openNormalBase K.absGal) (hNM : N ≤ M) (v : ℝ) :
    ((absGalStage K N v : Subgroup K.absGal) : Set K.absGal) * (M : Set K.absGal)
      = ((absGalStage K M v : Subgroup K.absGal) : Set K.absGal) := by
  rw [absGalStage_eq_stage (nonempty_stageGenerator K hN),
    absGalStage_eq_stage (nonempty_stageGenerator K hM)]
  exact stage_mul_coe_eq h K _ _ hNM v

/-- ★★★★★**`RamificationFiltration p`** —— 仮定は `UnramifiedBaseChangeHyp p` ただ 1 つ。

★★`Skeleton/PGC/Section2.lean` の `prop_2_1` / `prop_2_2` はこれを入力に取れる。 -/
noncomputable def ramificationFiltrationOfUnramifiedBaseChange (h : UnramifiedBaseChangeHyp p) :
    RamificationFiltration p :=
  ramificationFiltrationOfCompat (fun K => absGalStage_compat h K)

/-- ★★**構成した `Γ_K^v` は各有限段へ全射する** —— `Γ_K^v · M = Gal(L/L₀)^v` の引き戻し。
★これが「本物である」ことの証拠のひとつ(退化した族ではこの式は成り立たない)。 -/
theorem coe_ramificationFiltrationOfUnramifiedBaseChange_mul_coe (h : UnramifiedBaseChangeHyp p)
    (K : PAdicLocalField p) {M : Subgroup K.absGal} (hM : M ∈ openNormalBase K.absGal) (v : ℝ) :
    (((ramificationFiltrationOfUnramifiedBaseChange h).Gv K v : Subgroup K.absGal) :
          Set K.absGal) * (M : Set K.absGal)
      = ((absGalStage K M v : Subgroup K.absGal) : Set K.absGal) :=
  ramificationFiltrationOfStages_coe_mul_coe _ K hM v

/-! ### ★退化の自己検査 —— 近似 2 つと値が食い違うこと -/

/-- ★★**Y19b+c の `inertiaStageFiltration`(上からの近似)と一致しない**。 -/
theorem exists_stage_ne_inertiaStageFiltration (h : UnramifiedBaseChangeHyp p)
    (K : PAdicLocalField p) {N : Subgroup K.absGal} (hN : N ∈ openNormalBase K.absGal)
    (hne : ¬ absInertia K ≤ N) :
    ∃ v : ℝ, (absGalStageFiltration K (absGalStage_compat h K)).S N v
      ≠ (inertiaStageFiltration K).S N v := by
  simpa using exists_absGalStage_ne_inertiaStageFiltration K hN hne

/-- ★★**Y19 の `trivialStageFiltration`(下からの退化)とも一致しない**
(`v = 0` で `I_K ⊔ N` と `⊤` が違う)。 -/
theorem stage_ne_trivialStageFiltration (h : UnramifiedBaseChangeHyp p) (K : PAdicLocalField p)
    {N : Subgroup K.absGal} (hN : N ∈ openNormalBase K.absGal)
    (hne : absInertia K ⊔ N ≠ ⊤) :
    (absGalStageFiltration K (absGalStage_compat h K)).S N 0
      ≠ (trivialStageFiltration K.absGal).S N 0 := by
  rw [absGalStageFiltration_S, absGalStage_of_nonpos K hN le_rfl]
  show absInertia K ⊔ N ≠ (if (0 : ℝ) ≤ 0 then (⊤ : Subgroup K.absGal) else N)
  simpa using hne

/-! ## §9 `addVal` と `Ideal.ramificationIdx` の橋

★`QuotientRamIndexAverage.lean` が「mathlib の `Ideal.ramificationIdx` は `addVal` と
繋がっていない」と書いた穴をここで塞ぐ。★中身は `sSup {n | n ≤ v} = v` だけである。 -/

section AddValBridge

variable {A B : Type*} [CommRing A] [IsDomain A] [IsDiscreteValuationRing A]
  [CommRing B] [IsDomain B] [IsDiscreteValuationRing B]

/-- `addVal_map_eq_mul_addVal` を `Algebra` 無しの環準同型で使う形。 -/
theorem addVal_ringHom_mul (φ : A →+* B) {ϖ : A} (hϖ : Irreducible ϖ)
    (he : addVal B (φ ϖ) ≠ 0) (a : A) : addVal B (φ a) = addVal B (φ ϖ) * addVal A a := by
  letI := φ.toAlgebra
  show addVal B (algebraMap A B a) = addVal B (algebraMap A B ϖ) * addVal A a
  exact addVal_map_eq_mul_addVal hϖ he a

/-- `addVal a = 1` なら `a` は素元(`addVal_uniformizer` の逆)。 -/
theorem irreducible_of_addVal_eq_one {a : A} (h : addVal A a = 1) : Irreducible a := by
  obtain ⟨ϖ, hϖ⟩ := IsDiscreteValuationRing.exists_irreducible A
  have ha : a ≠ 0 := by
    intro h0
    rw [h0, addVal_zero] at h
    exact (by decide : (⊤ : ℕ∞) ≠ 1) h
  obtain ⟨n, u, rfl⟩ := IsDiscreteValuationRing.eq_unit_mul_pow_irreducible ha hϖ
  rw [addVal_def' u hϖ n] at h
  have hn : n = 1 := by exact_mod_cast h
  subst hn
  exact Associated.irreducible ⟨u, mul_comm _ _⟩ (by simpa using hϖ)

end AddValBridge

theorem injective_adjoinIntegersRingHom (K : PAdicLocalField p) {x x' : K.closure}
    (hle : IntermediateField.adjoin K.carrier ({x} : Set K.closure)
      ≤ IntermediateField.adjoin K.carrier ({x'} : Set K.closure)) :
    Function.Injective (adjoinIntegersRingHom K hle) := by
  intro z₁ z₂ hz
  refine Subtype.ext (Subtype.ext ?_)
  have h1 := congrArg
    (fun w : adjoinIntegers K x' =>
      ((w : IntermediateField.adjoin K.carrier ({x'} : Set K.closure)) : K.closure)) hz
  simpa only [adjoinIntegersRingHom_apply, coe_adjoinIntegersIncl] using h1

theorem val_algebraMap_adjoinIntegers (K : PAdicLocalField p) (x : K.closure)
    (y : 𝒪[K.carrier]) :
    (((algebraMap 𝒪[K.carrier] (adjoinIntegers K x) y : adjoinIntegers K x) :
        IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure)
      = algebraMap K.carrier K.closure (y : K.carrier) := rfl

/-- ★★**`addVal` と `Ideal.ramificationIdx` の橋** —— `𝒪_K` の素元の `𝒪_{K(x)}` での
付値はちょうど分岐指数 `e(K(x)/K)` である。

★`Ideal.ramificationIdx p P = sSup {n | map p ≤ P^n}` の定義に戻って、その集合が
`Set.Iic v` であることを見るだけ(`csSup_Iic`)。 -/
theorem addVal_algebraMap_eq_ramificationIndex (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    {ϖ : 𝒪[K.carrier]} (hϖ : Irreducible ϖ) :
    addVal (adjoinIntegers K x) (algebraMap 𝒪[K.carrier] (adjoinIntegers K x) ϖ)
      = (ramificationIndex K x : ℕ∞) := by
  haveI := isDiscreteValuationRing_carrierIntegers K
  obtain ⟨π, hπ⟩ := IsDiscreteValuationRing.exists_irreducible (adjoinIntegers K x)
  have hne : algebraMap 𝒪[K.carrier] (adjoinIntegers K x) ϖ ≠ 0 := by
    intro h0
    exact hϖ.ne_zero (injective_algebraMap_adjoinIntegers K x (by rw [h0, map_zero]))
  obtain ⟨v, hv⟩ : ∃ v : ℕ,
      addVal (adjoinIntegers K x) (algebraMap 𝒪[K.carrier] (adjoinIntegers K x) ϖ)
        = (v : ℕ∞) := by
    obtain ⟨v, hv⟩ := ENat.ne_top_iff_exists.mp (fun h => hne (addVal_eq_top_iff.1 h))
    exact ⟨v, hv.symm⟩
  have hset : {n : ℕ | Ideal.map (algebraMap 𝒪[K.carrier] (adjoinIntegers K x))
      (maximalIdeal 𝒪[K.carrier]) ≤ (maximalIdeal (adjoinIntegers K x)) ^ n} = Set.Iic v := by
    ext n
    rw [(IsDiscreteValuationRing.irreducible_iff_uniformizer ϖ).1 hϖ, Ideal.map_span]
    simp only [Set.image_singleton, Set.mem_setOf_eq, Ideal.span_le, Set.singleton_subset_iff,
      SetLike.mem_coe, Set.mem_Iic]
    rw [mem_maximalIdeal_pow_iff_le_addVal hπ n, hv]
    exact Nat.cast_le
  rw [hv, ramificationIndex, Ideal.ramificationIdx, hset, csSup_Iic]

/-! ## §10 ★★★`e = |I|` —— 分岐指数は惰性群の位数

★★これが `compat` に残っていた**唯一の数学**である。

段取り(すべて既存在庫の組み合わせ):
* `Nat.card I(L/K) = [L : L₀]`(`IsGalois.card_fixingSubgroup_eq_finrank` +
  Y19b+c の `inertiaGal_eq_fixingSubgroup`)。
* `[L:K] = [L₀:K] · [L:L₀]`(`Module.finrank_mul_finrank`)。
* ★**`[L₀:K] = f(L/K)`** —— `L₀ = K^ur ⊓ L`(Y19b+c `lift_fixedField_inertiaGal`)を
  原始元で `K(w)` と書き、`IsUnramifiedAdjoin K w` から `[K(w):K] = f(K(w)/K) ≤ f(L/K)`、
  逆向きは `AbelianSplitUnramified.lean` の `exists_unramified_subextension`
  (最大不分岐部分拡大の次数が `f`)。
* `e·f = [L:K]`(`ramificationIndex_mul_inertiaDegree`)で `f` を約す。 -/

/-- `L/K` は正規かつ(標数 0 なので)分離的、すなわち Galois。 -/
theorem isGalois_adjoin (K : PAdicLocalField p) (x : K.closure)
    [Normal K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))] :
    IsGalois K.carrier ↥(IntermediateField.adjoin K.carrier ({x} : Set K.closure)) := by
  haveI := isGalois_closure K
  haveI : Algebra.IsSeparable K.carrier
      ↥(IntermediateField.adjoin K.carrier ({x} : Set K.closure)) :=
    Algebra.isSeparable_tower_bot_of_isSeparable K.carrier
      ↥(IntermediateField.adjoin K.carrier ({x} : Set K.closure)) K.closure
  exact ⟨⟩

set_option maxHeartbeats 1000000 in
/-- ★★★**`[L₀ : K] = f(L/K)`** —— 惰性群の固定体はちょうど最大不分岐部分拡大である。 -/
theorem finrank_fixedField_inertiaGal (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    [Normal K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))] :
    Module.finrank K.carrier ↥(IntermediateField.fixedField (inertiaGalAdjoin K x))
      = inertiaDegree K x := by
  haveI := isGalois_closure K
  set L := IntermediateField.adjoin K.carrier ({x} : Set K.closure) with hLdef
  set F := IntermediateField.fixedField (inertiaGalAdjoin K x) with hFdef
  set M := IntermediateField.lift F with hMdef
  have hM : M = unramifiedClosure K ⊓ L := lift_fixedField_inertiaGal K L
  have hML : M ≤ L := by rw [hM]; exact inf_le_right
  have hMur : M ≤ unramifiedClosure K := by rw [hM]; exact inf_le_left
  haveI : FiniteDimensional K.carrier ↥F := IntermediateField.finiteDimensional_left F
  have hequiv : ↥F ≃ₗ[K.carrier] ↥M := (IntermediateField.equivMap F L.val).toLinearEquiv
  haveI : FiniteDimensional K.carrier ↥M := hequiv.finiteDimensional
  have hFM : Module.finrank K.carrier ↥F = Module.finrank K.carrier ↥M := hequiv.finrank_eq
  haveI : Algebra.IsSeparable K.carrier ↥M :=
    Algebra.isSeparable_tower_bot_of_isSeparable K.carrier ↥M K.closure
  obtain ⟨w, hw⟩ : ∃ w : K.closure,
      IntermediateField.adjoin K.carrier ({w} : Set K.closure) = M := by
    obtain ⟨α, hα⟩ := Field.exists_primitive_element K.carrier ↥M
    refine ⟨(α : K.closure), ?_⟩
    have h1 := IntermediateField.adjoin_map K.carrier ({α} : Set ↥M) M.val
    rw [hα] at h1
    have h2 : IntermediateField.map M.val (⊤ : IntermediateField K.carrier ↥M) = M :=
      IntermediateField.lift_top K.carrier M
    rw [h2] at h1
    have h3 : IntermediateField.adjoin K.carrier ({(α : K.closure)} : Set K.closure)
        = IntermediateField.adjoin K.carrier (⇑M.val '' {α}) := by
      congr 1
      simp
    exact h3.trans h1.symm
  have hu : IsUnramifiedAdjoin K w :=
    (adjoin_le_unramifiedClosure_iff K _).1 (by rw [hw]; exact hMur)
  have hup : Module.finrank K.carrier ↥M ≤ inertiaDegree K x := by
    have h3 : Module.finrank K.carrier ↥M = inertiaDegree K w := by
      rw [← hw]
      exact (inertiaDegree_eq_finrank_of_isUnramified K _ hu).symm
    rw [h3]
    exact inertiaDegree_le_of_adjoin_le K (by rw [hw]; exact hML)
  obtain ⟨z, huz, hfz, hzx⟩ := exists_unramified_subextension K x
  have hzM : IntermediateField.adjoin K.carrier ({z} : Set K.closure) ≤ M := by
    rw [hM]
    exact le_inf (adjoin_le_unramifiedClosure K huz) hzx
  have hlow : inertiaDegree K x ≤ Module.finrank K.carrier ↥M := by
    rw [← hfz]
    exact Submodule.finrank_mono (R := K.carrier) (M := K.closure)
      (s := (IntermediateField.adjoin K.carrier ({z} : Set K.closure)).toSubalgebra.toSubmodule)
      (t := M.toSubalgebra.toSubmodule) (fun y hy => hzM hy)
  rw [hFM]
  exact le_antisymm hup hlow

set_option maxHeartbeats 1000000 in
/-- ★★★★**`e(K(x)/K) = |I(K(x)/K)|`**。 -/
theorem ramificationIndex_eq_card_inertiaGal (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    [Normal K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))] :
    ramificationIndex K x = Nat.card ↥(inertiaGalAdjoin K x) := by
  haveI := isGalois_adjoin K x
  have h1 : Nat.card ↥(inertiaGalAdjoin K x)
      = Module.finrank ↥(IntermediateField.fixedField (inertiaGalAdjoin K x))
          ↥(IntermediateField.adjoin K.carrier ({x} : Set K.closure)) := by
    have h := IsGalois.card_fixingSubgroup_eq_finrank
      (IntermediateField.fixedField (inertiaGalAdjoin K x))
    rw [← inertiaGal_eq_fixingSubgroup] at h
    exact h
  have h2 := Module.finrank_mul_finrank K.carrier
    ↥(IntermediateField.fixedField (inertiaGalAdjoin K x))
    ↥(IntermediateField.adjoin K.carrier ({x} : Set K.closure))
  rw [finrank_fixedField_inertiaGal K x, ← h1,
    ← ramificationIndex_mul_inertiaDegree K x] at h2
  have hfpos : 0 < inertiaDegree K x := Nat.pos_of_ne_zero (fun h0 => by
    have hpos : 0 < Module.finrank K.carrier
        ↥(IntermediateField.adjoin K.carrier ({x} : Set K.closure)) := Module.finrank_pos
    rw [← ramificationIndex_mul_inertiaDegree K x, h0, mul_zero] at hpos
    exact lt_irrefl 0 hpos)
  rw [mul_comm (ramificationIndex K x) (inertiaDegree K x)] at h2
  exact (Nat.eq_of_mul_eq_mul_left hfpos h2).symm

/-- ★★★**絶対分岐指数は惰性群の位数**(§9 の橋 + §10)。 -/
theorem addVal_algebraMap_eq_card_inertiaGal (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    [Normal K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    {ϖ : 𝒪[K.carrier]} (hϖ : Irreducible ϖ) :
    addVal (adjoinIntegers K x) (algebraMap 𝒪[K.carrier] (adjoinIntegers K x) ϖ)
      = (Nat.card ↥(inertiaGalAdjoin K x) : ℕ∞) := by
  rw [addVal_algebraMap_eq_ramificationIndex K x hϖ, ramificationIndex_eq_card_inertiaGal K x]

/-! ## §11 ★★★★★不分岐底変換の証明と `RamificationFiltration p` -/

set_option maxHeartbeats 400000 in
/-- ★★★★★**不分岐底変換で素元は素元のまま** —— `compat` に残っていた仮定が定理になった。

段取り(`E := v_{𝒪_{L′}}(π)`、`H := ker(I(L′/K) ↠ I(L/K))`):
* `v_{𝒪_{L′}}(ι a) = E · v_{𝒪_L}(a)`(§9 `addVal_ringHom_mul`)。
* `𝒪_K` の素元 `ϖ_K` に §10 を 2 回当てて `|I(L′/K)| = E · |I(L/K)|`、
  `|I(L′/K)| = |H| · |I(L/K)|`(§1 `card_eq_card_ker_mul`)と比べて **`E = |H|`**。
* `v_{𝒪_{L′}}(ι_C c) = |H| · v_C(c)`(`FixedRingTower.lean`
  `addVal_map_eq_card_mul_fixedRing`)に `c := stageBaseHom π` を入れて `v_C(c) = 1`。
* `addVal = 1 → Irreducible`(§9)。 -/
theorem unramifiedBaseChange : UnramifiedBaseChangeHyp p := by
  intro K x x' _ _ _ _ hle
  haveI : Fintype ↥(inertiaRestrict K hle).ker := Fintype.ofFinite _
  haveI := isDiscreteValuationRing_carrierIntegers K
  have hπ : Irreducible (stageUniformizer K x) := irreducible_stageUniformizer K
  have hne : adjoinIntegersRingHom K hle (stageUniformizer K x) ≠ 0 := by
    intro h0
    exact hπ.ne_zero (injective_adjoinIntegersRingHom K hle (by rw [h0, map_zero]))
  have hE0 : addVal (adjoinIntegers K x')
      (adjoinIntegersRingHom K hle (stageUniformizer K x)) ≠ 0 := by
    intro h0
    exact hπ.not_isUnit (isUnit_of_map_unit (adjoinIntegersRingHom K hle) _
      (addVal_eq_zero_iff.1 h0))
  have hEtop : addVal (adjoinIntegers K x')
      (adjoinIntegersRingHom K hle (stageUniformizer K x)) ≠ ⊤ :=
    fun h0 => hne (addVal_eq_top_iff.1 h0)
  obtain ⟨e, he⟩ : ∃ e : ℕ, addVal (adjoinIntegers K x')
      (adjoinIntegersRingHom K hle (stageUniformizer K x)) = (e : ℕ∞) := by
    obtain ⟨e, he⟩ := ENat.ne_top_iff_exists.mp hEtop
    exact ⟨e, he.symm⟩
  obtain ⟨ϖK, hϖK⟩ := IsDiscreteValuationRing.exists_irreducible 𝒪[K.carrier]
  have hcompat : algebraMap 𝒪[K.carrier] (adjoinIntegers K x') ϖK
      = adjoinIntegersRingHom K hle (algebraMap 𝒪[K.carrier] (adjoinIntegers K x) ϖK) := by
    refine Subtype.ext (Subtype.ext ?_)
    simp only [adjoinIntegersRingHom_apply, coe_adjoinIntegersIncl,
      val_algebraMap_adjoinIntegers]
  have hcard : Nat.card ↥(inertiaGalAdjoin K x')
      = Nat.card ↥(inertiaRestrict K hle).ker * Nat.card ↥(inertiaGalAdjoin K x) :=
    card_eq_card_ker_mul (surjective_inertiaRestrict K hle)
  have h3 : (Nat.card ↥(inertiaGalAdjoin K x') : ℕ∞)
      = (e : ℕ∞) * (Nat.card ↥(inertiaGalAdjoin K x) : ℕ∞) := by
    rw [← addVal_algebraMap_eq_card_inertiaGal K x' hϖK, hcompat,
      addVal_ringHom_mul _ hπ hE0, he, addVal_algebraMap_eq_card_inertiaGal K x hϖK]
  have hE : e = Nat.card ↥(inertiaRestrict K hle).ker := by
    have h4 : Nat.card ↥(inertiaGalAdjoin K x')
        = e * Nat.card ↥(inertiaGalAdjoin K x) := by exact_mod_cast h3
    rw [hcard] at h4
    exact (Nat.eq_of_mul_eq_mul_right Nat.card_pos h4.symm)
  have h7 := addVal_map_eq_card_mul_fixedRing
    (A := ↥(fixedRing (adjoinIntegers K x') (inertiaGalAdjoin K x')))
    (H := (inertiaRestrict K hle).ker)
    (irreducible_stageUniformizer K)
    (fun b => exists_sub_mem_maximalIdeal_inertiaFixedRing K x' b)
    (adjoin_stageUniformizer_eq_top K x')
    (stageBaseHom K hle (stageUniformizer K x))
  have hbase : algebraMap ↥(fixedRing (adjoinIntegers K x') (inertiaRestrict K hle).ker)
        (adjoinIntegers K x') (stageBaseHom K hle (stageUniformizer K x))
      = adjoinIntegersRingHom K hle (stageUniformizer K x) :=
    coe_stageBaseHom K hle (stageUniformizer K x)
  rw [hbase, he, hE] at h7
  have hHpos : 0 < Nat.card ↥(inertiaRestrict K hle).ker := Nat.card_pos
  have hvtop : addVal ↥(fixedRing (adjoinIntegers K x') (inertiaRestrict K hle).ker)
      (stageBaseHom K hle (stageUniformizer K x)) ≠ ⊤ := by
    intro h0
    rw [h0, ENat.mul_top (by exact_mod_cast hHpos.ne')] at h7
    exact (ENat.coe_ne_top _) h7
  obtain ⟨m, hm⟩ : ∃ m : ℕ, addVal ↥(fixedRing (adjoinIntegers K x') (inertiaRestrict K hle).ker)
      (stageBaseHom K hle (stageUniformizer K x)) = (m : ℕ∞) := by
    obtain ⟨m, hm⟩ := ENat.ne_top_iff_exists.mp hvtop
    exact ⟨m, hm.symm⟩
  rw [hm, ← Nat.cast_mul] at h7
  have h8 : Nat.card ↥(inertiaRestrict K hle).ker * 1
      = Nat.card ↥(inertiaRestrict K hle).ker * m := by
    rw [mul_one]; exact_mod_cast h7
  have hm1 : m = 1 := (Nat.eq_of_mul_eq_mul_left hHpos h8).symm
  exact irreducible_of_addVal_eq_one (by rw [hm, hm1, Nat.cast_one])

/-- ★★★★★★**`Interface/PGC/LocalFieldData.lean` の `RamificationFiltration p` の本物**。

★★**仮定は 1 つも無い。** `Skeleton/PGC/Section2.lean` の `prop_2_1` / `prop_2_2` は
これをそのまま入力に取れる。 -/
noncomputable def ramificationFiltration (p : ℕ) [Fact p.Prime] : RamificationFiltration p :=
  ramificationFiltrationOfUnramifiedBaseChange unramifiedBaseChange

/-- ★★**構成した `Γ_K^v` は各有限段へ全射する**(仮定なしの形)。 -/
theorem coe_ramificationFiltration_mul_coe (K : PAdicLocalField p) {M : Subgroup K.absGal}
    (hM : M ∈ openNormalBase K.absGal) (v : ℝ) :
    (((ramificationFiltration p).Gv K v : Subgroup K.absGal) : Set K.absGal)
        * (M : Set K.absGal)
      = ((absGalStage K M v : Subgroup K.absGal) : Set K.absGal) :=
  coe_ramificationFiltrationOfUnramifiedBaseChange_mul_coe unramifiedBaseChange K hM v

end ABC3.Found.PGC
