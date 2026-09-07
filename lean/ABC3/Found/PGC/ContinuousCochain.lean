import ABC3.Found.PGC.CohomologyColimit

/-!
# [pGC] Proposition 1.1 —— 連続コチェイン複体 `C^n_cont(G,A)` と `H^n_cont(G,A)`

`Found/PGC/CohomologyColimit.lean` は塔 `S i ↓ 1` に沿った colimit 側
(`H2Colim` / 比較射 `H2ColimToH2` / その単射性)を離散側で完全に立てた。
★そこで残ったのが**相手側の対象が無い**という一点である:

> mathlib には**連続コホモロジーが次数 0 以外に無い**((c)-2 の測定、本ファイルで再測)。
> したがって `colim_i H^n(G ⧸ S i, A^{S i}) ≅ H^n_cont(G,A)` を述べる**相手が存在しない**。

本ファイルはその相手 —— **連続コチェイン複体 `C^n_cont(G,A)` と `H^n_cont(G,A)`** ——
を作る。★次数 1 の雛形 `smoothCocycles₁`(`CohomologyColimit.lean`)を次数 `n` に上げたものである。

## ★★★★冒頭に明示する —— 離散と連続の違いは結論を変える

★★`CohomologyColimit.lean` が実証したとおり、**離散側では `colim_i H¹(S i,A) = 0` は偽**である
(反例 `G = Ẑ`, `A = ℚ/ℤ`)。連続でない 1-コサイクルはどの `nẐ` の上でも消えない。
★したがって「連続」は飾りではなく、**主張の真偽を決める**。本ファイルの `C^n_cont` は
その「連続」を Lean の中で初めて次数 `n` で定義したものである。

## ★★★本ファイルの到達点

1. ★★★**`C^n_cont(G,A)`**(`contCochains`)——
   ★**位相を使わない述語**(`LocConstNearOne`)で定義した。
   `𝓝 (1 : G)` しか使わず、係数 `A` の位相は**一切要求しない**
   (★`VanishesNearOne` の流儀そのまま)。
2. ★★★**`H^n_cont(G,A)`**(`Hcont`)—— `Z^n_cont ⧸ B^n_cont`。
   ★`d` が `C^n_cont` を保つこと(`d_mem_contCochains`)が本ファイルの中核で、
   ★中身は 3 本の抽象核(`LocConstNearOne.comp` / `isNearOneEquivariant_contractNth` /
   `LocConstNearOne.rhoTerm`)に入っている。
3. ★★★(c)-3 の 3 番目(colimit との比較)は**コチェイン水準まで**入った:
   `inflCochain`(inflation)、★`inflCochain_mem_contCochains`(**`S` が `1` の近傍なら
   inflation は連続コチェインを与える**)、`inflHcont`(比較射)、
   ★`inflHcont_inflStepCochain`(**塔と両立する**)。
   ★止まった場所は下の「★★(c)-4 との境界」に明記した。

## ★★★退化していない証拠(この 4 本で確かめた)

* ★★★**`LocConstNearOne f ↔ Continuous f`**(`locConstNearOne_iff_continuous`、
  係数が離散・添字が有限)。★**両向き**を証明した。したがって `C^n_cont` は
  「連続コチェイン」そのものであり、痩せても太ってもいない
  (`mem_contCochains_iff_continuous`)。
* ★★★**次数 1 に潰すと `smoothCocycles₁` に一致する**
  (`mem_contCochains_one_iff_mem_smoothCocycles₁`)。★1-コサイクルについては
  「なめらか(本ファイル)」と「`1` の近傍で消える(前ファイル)」は**同値**である。
  ★★仮定は**コサイクル条件だけ**で、副有限性も離散性も要らない。
* ★★**`C^0_cont(G,A) = A`、`Z^0_cont(G,A) = A^G`**(`contCochains_zero` /
  `mem_contCocycles_zero_iff`)。したがって `H^0_cont(G,A) = A^G`(`HcontZeroEquiv`)。
  ★「なめらかベクトル」に痩せていない(次数 0 に条件が付かない)ことの確認である。
* ★**`G` が離散なら `C^n_cont = C^n`**(`contCochains_of_discreteTopology`)。
* ★次数 1 のコサイクル条件が mathlib の `cocycles₁` と一致する(`mem_ker_d_one_iff`)。

## ★★★`A` の離散性をどこで使ったか(★5 本だけ)

★`DiscreteTopology` を使うのは次の 5 本**だけ**である。
`continuous_of_locConstNearOne` / `locConstNearOne_of_continuous` /
`locConstNearOne_iff_continuous` / `smoothAction_iff_continuousAt` と、
その具体層 `mem_contCochains_iff_continuous` / `smoothAction_repAddHom_iff`。
★★**`C^n_cont` と `H^n_cont` の定義にも、`d` が保たれることの証明にも、離散性は要らない。**
★これは前ファイルで `vanishesNearOne_of_continuousAt` 1 本だけが離散性を使ったのと同じ設計である。

## ★★設計 —— 抽象核と具体層

★★**§1〜§5 は「コホモロジー・分岐・付値・Galois の語彙が 1 語も出ない核」**である。
出てくるのは群 `G`、その位相、加法群 `M`、`ρ : G → M →+ M`(#209)だけ。

* **§1** `LocConstNearOne` —— 「各点で右移動について局所定数」。`smoothFuns`(部分加群)。
* **§2** `IsNearOneEquivariant` と `LocConstNearOne.comp` ——
  ★★**座標変換で保たれる**ための十分条件。
  ★`isNearOneEquivariant_contractNth`(隣り合う 2 座標の積)が唯一の非自明な核で、
  ★使うのは**位相群での共役の連続性**だけ(`exists_nhds_one_split` + `ContinuousAt`)。
* **§3** `SmoothAction` と `LocConstNearOne.rhoTerm` —— 微分の「作用の項」。
* **§4** 次数 1 の核(`locConstNearOne_one_iff_vanishes`)——
  ★**1-コサイクルなら「なめらか」= 「`1` の近傍で消える」**。
* **§5** 離散性の核 —— `LocConstNearOne` ↔ `Continuous`。

具体層(§6 以降)は核に代入するだけである。★実測: 核は 0.09〜0.37 秒、具体層は 0.18〜0.77 秒。

## ★★★測定 —— mathlib に何があり何が無いか(2026-09-07、コマンドと出力)

```
awk -F'\t' '$2 ~ /inhomogeneousCochains/ {print $2"\t"$3}' .cache/mathlib-index.txt
```
→ 6 件。`inhomogeneousCochains.d`(`Basic.lean:86`、★**根の名前空間**)、
`groupCohomology.inhomogeneousCochains.d_comp_d`(`:142`)、
`groupCohomology.inhomogeneousCochains.d_def`(`:137`)ほか。
★**次数 `n` の微分と `d ∘ d = 0` は在る**(本ファイルはこれを使う。自分で書いていない)。

```
awk -F'\t' '$2 ~ /^groupCohomology\.(cochainsMap|cocyclesMap|π|cocycles)/ ...' .cache/mathlib-index.txt
```
→ `cochainsMap`(`Functoriality.lean:52`)、`cocycles`(`Basic.lean:159`)、
`cocyclesMk`(`:163`)、`π`(`:192`)、`iCocycles`。
★**`cochainsMap` が在るので inflation の可換性(`inflCochain_comm_d`)は 3 行で済んだ。**

```
awk -F'\t' '$2 ~ /^(Fin\.contractNth|exists_nhds_one_split|Filter\.iInter_mem|...)/ ...'
```
→ `Fin.contractNth_apply_of_lt/of_eq/of_gt`(`Data/Fin/Tuple/Basic.lean:1254-1263`)、
`exists_nhds_one_split`(`Topology/Algebra/Monoid.lean:581`)、
`Filter.iInter_mem`(`Order/Filter/Finite.lean:48`)、
`Submodule.quotEquivOfEqBot`(`LinearAlgebra/Quotient/Basic.lean:399`)。★全部在った。

★★**無いもの(本ファイルが作った理由)**: **次数 1 以上の連続コホモロジー**。
(c)-2 の測定と同じ結論で、本ファイルはそこを埋める。

★重複検査(#158): 本ファイルの 60 個の名前を実名前空間 `ABC3.Found.PGC` に
`addToEnv` で積み、`already declared` が出るか見た。
→ ★**1 件**(`repAddHom_mul`。`Transgression.lean` に同じ主張が既にあったので**自分の版を捨てた**)。
残り 59 件は衝突なし。さらに `.cache/decl-index.txt` と `.cache/mathlib-index.txt` を
19 個の主要名で grep して 0 件。

## ★正規性・可判定性・選択公理をどこで使ったか

★**正規性を使う宣言**: §10 の inflation 系だけ(`G ⧸ S` を作るため)——
`inflCochain` / `inflCochain_apply` / `inflCochain_comm_d` /
`inflCochain_mem_contCochains` / `inflCochain_mem_contCocycles` /
`inflCochain_d_mem_contCoboundaries` / `inflStepCochain` /
`inflStepCochain_comm_d` / `inflCochain_inflStepCochain` /
`inflHcont` / `inflHcont_d_eq_zero` / `inflHcont_inflStepCochain`。
★★**`C^n_cont` と `H^n_cont` の定義・性質は正規性を一切使わない。**

★**可判定性**: 本ファイルは `Module.DirectLimit` を作らないので `[DecidableEq ι]` は
**1 つも要らない**(#221)。`Classical.decEq` も差し込んでいない。

★**選択公理**: `#print axioms` は
`LocConstNearOne` の核 4 本(`.add` / `.smul` / `.comp` / `.rhoTerm`)と
`isNearOneEquivariant_contractNth` / `locConstNearOne_one_iff` /
`locConstNearOne_one_iff_vanishes` が ★**`[propext, Quot.sound]` のみ**。
`Rep` や `ModuleCat` を触るもの(§6 以降)は `[propext, Classical.choice, Quot.sound]`。

## ★★(c)-4 との境界 —— どこで止まったか

★**入ったのはここまで**: `inflHcont : Z^n(G ⧸ S, A^S) →ₗ H^n_cont(G,A)` と
「コバウンダリを殺す」(`inflHcont_d_eq_zero`)と「塔と両立する」
(`inflHcont_inflStepCochain`)。★**比較射の中身はすべて揃っている。**

★**止まった理由は 1 つだけ**である。

> mathlib の `groupCohomology B n`(`n ≥ 3`)は**ホモロジー対象**であって
> `ker ⧸ im` の形では与えられていない。`groupCohomology.π B n : cocycles B n ⟶ H^n` は
> 在るが、**そこから出る射を作る descent(余核の普遍性)** を経由しないと
> `H^n(G ⧸ S, A^S) → H^n_cont(G,A)` にはならない。
> ★次数 2 なら `H2π` と `H2_induction_on` が在る(`CohomologyColimit.lean` が使っている)ので、
> **次数 2 に限れば `H2Colim → H^2_cont` は本ファイルの在庫だけで作れる**と見込む。
> ★そのとき必要なのは `cocycles₂ B`(`G × G → B`)と `ker (d B 2)`(`(Fin 2 → G) → B`)の
> 突き合わせで、本ファイルの `mem_ker_d_one_iff` の次数 2 版(★`Fin 3` の 3 つの
> `contractNth` を計算する)である。★次数 1 版は 30 行で済んだ。

★**(c)-4 の見通し**: `H2ColimToH2` の**全射性**は本ファイルでは動かない。
★ただし本ファイルの `inflCochain_mem_contCochains` は
「**`H²(G,A)` ではなく `H²_cont(G,A)` を相手にすれば全射性が期待できる**」ことを示している
(連続 2-コサイクルは開正規部分群を経由するから)。★相手を差し替える必要がある。

## 逸脱の記録

* **逸脱 1**: 原典 (pGC) は `Γ_K` の連続コホモロジーを使う。本ファイルはそれを
  **mathlib の非斉次コチェイン複体の部分複体**として定義した。
  ★`Transgression.lean` / `CohomologyColimit.lean` の逸脱 1(離散で書く)を、
  ★**本ファイルで初めて連続側に戻した**。
* **逸脱 2**: 「連続」を位相ではなく述語 `LocConstNearOne` で書いた。
  ★係数の位相を要求しないためである。★同値性(`locConstNearOne_iff_continuous`)を
  **両向き**証明したので、逸脱は名前だけである。
* **逸脱 3**: すべての宣言で `Rep.{0} k G` と宇宙を固定した。
  ★`groupCohomology.inhomogeneousCochains.d_comp_d` が `Rep.{u,u,u}` を要求する一方
  `inhomogeneousCochains.d` は係数の宇宙が自由なため、`Rep k G` のままだと
  `Application type mismatch: A has type Rep.{u_1, 0, 0} k G but is expected to have type
  Rep.{u_1, u_1, u_1}` になる(#224)。★数学的な制限ではない。
* **逸脱 4**: 塔の段 `S` に課すのは `(S : Set G) ∈ 𝓝 1` だけで、開性も正規性も
  必要な所でしか使わない。★原典の `Γ_K` の開正規部分群はこれを満たす
  (`profiniteTower_mem_nhds_one`)。
-/

namespace ABC3.Found.PGC

open CategoryTheory groupCohomology

/-- 原典の位置(★`CohomologyColimit.lean` の `cohomologyColimit.src` と同じ項目)。

原文 (pGC p.3):
> Proposition 1.1: The cyclotomic character χ : ΓK → Zp× can be recovered entirely
> group-theoretically from ΓK.

★原典は連続コホモロジーを「well-known」で畳んでいる。本ファイルは mathlib 側の欠落
(連続コホモロジーが次数 0 しか無い)を、次数 `n` の連続コチェイン複体で埋める位置づけである。 -/
def continuousCochain.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 3, item := "Proposition 1.1", sectionId := "prop-1-1" }

/-! ## 1. ★★★抽象核 I —— なめらかな関数(語彙ゼロ)

★出てくるのは群 `G`、その位相、加法群 `M` だけ。
★コホモロジー・分岐・付値・Galois の語彙は 1 語も無い。 -/

section AbstractCoreSmoothFun

open Topology

variable {G : Type*} [Group G] [TopologicalSpace G] {M : Type*} [AddCommGroup M]

/-- ★★★核 —— `f : (ι → G) → M` が**各点で右移動について局所定数**であること。

★これが本ファイルの「連続コチェイン」の定義である。★`M` の位相は要求しない。
★`M` が離散で `ι` が有限なら `Continuous f` と同値(`locConstNearOne_iff_continuous`)。 -/
def LocConstNearOne {ι : Type*} (f : (ι → G) → M) : Prop :=
  ∀ g : ι → G, ∃ V ∈ 𝓝 (1 : G), ∀ h : ι → G, (∀ i, h i ∈ V) → f (fun i => g i * h i) = f g

omit [AddCommGroup M] in
theorem locConstNearOne_const {ι : Type*} (a : M) :
    LocConstNearOne (fun _ : ι → G => a) :=
  fun _ => ⟨Set.univ, Filter.univ_mem, fun _ _ => rfl⟩

theorem LocConstNearOne.add {ι : Type*} {f₁ f₂ : (ι → G) → M}
    (h₁ : LocConstNearOne f₁) (h₂ : LocConstNearOne f₂) : LocConstNearOne (f₁ + f₂) := by
  intro g
  obtain ⟨V₁, hV₁, e₁⟩ := h₁ g
  obtain ⟨V₂, hV₂, e₂⟩ := h₂ g
  refine ⟨V₁ ∩ V₂, Filter.inter_mem hV₁ hV₂, fun h hh => ?_⟩
  show f₁ _ + f₂ _ = f₁ g + f₂ g
  rw [e₁ h fun i => (hh i).1, e₂ h fun i => (hh i).2]

theorem LocConstNearOne.smul {k : Type*} [Semiring k] [Module k M] {ι : Type*} (a : k)
    {f : (ι → G) → M} (hf : LocConstNearOne f) : LocConstNearOne (a • f) := by
  intro g
  obtain ⟨V, hV, e⟩ := hf g
  refine ⟨V, hV, fun h hh => ?_⟩
  show a • f _ = a • f g
  rw [e h hh]

/-- ★**なめらかな関数のなす部分加群**(核)。 -/
def smoothFuns (k : Type*) [Semiring k] [Module k M] (ι : Type*) :
    Submodule k ((ι → G) → M) where
  carrier := {f | LocConstNearOne f}
  add_mem' h₁ h₂ := h₁.add h₂
  zero_mem' := locConstNearOne_const 0
  smul_mem' a _ hf := hf.smul a

@[simp] theorem mem_smoothFuns {k : Type*} [Semiring k] [Module k M] {ι : Type*}
    {f : (ι → G) → M} : f ∈ smoothFuns (M := M) k ι ↔ LocConstNearOne f := Iff.rfl

end AbstractCoreSmoothFun

/-! ## 2. ★★★抽象核 II —— 座標変換で保たれること(語彙ゼロ)

★`d` の各項は「座標変換との合成」である。それが `LocConstNearOne` を保つ条件を切り出す。 -/

section AbstractCoreComp

open Topology

variable {G : Type*} [Group G] [TopologicalSpace G] {M : Type*} [AddCommGroup M]
  {ι κ : Type*}

/-- 核 —— 座標変換 `Φ` が「`1` の近くの右移動を `1` の近くの右移動に写す」。 -/
def IsNearOneEquivariant (Φ : (κ → G) → (ι → G)) : Prop :=
  ∀ (g : κ → G), ∀ V ∈ 𝓝 (1 : G), ∃ W ∈ 𝓝 (1 : G), ∀ h : κ → G, (∀ i, h i ∈ W) →
    ∀ i, (Φ g i)⁻¹ * Φ (fun i => g i * h i) i ∈ V

omit [AddCommGroup M] in
/-- ★★核 —— なめらかさは `IsNearOneEquivariant` な座標変換で保たれる。

★★「移動の差 `(Φ g i)⁻¹ * Φ (g·h) i` が小さい」とだけ言えばよい形にしたので、
座標変換の等式は `mul_inv_cancel_left` で**自動的に**出る(★原典の段取りより短い)。 -/
theorem LocConstNearOne.comp {Φ : (κ → G) → (ι → G)} (hΦ : IsNearOneEquivariant Φ)
    {f : (ι → G) → M} (hf : LocConstNearOne f) : LocConstNearOne (fun g => f (Φ g)) := by
  intro g
  obtain ⟨V, hV, e⟩ := hf (Φ g)
  obtain ⟨W, hW, hWV⟩ := hΦ g V hV
  refine ⟨W, hW, fun h hh => ?_⟩
  have key : Φ (fun i => g i * h i) = fun i => Φ g i * ((Φ g i)⁻¹ * Φ (fun i => g i * h i) i) := by
    funext i; rw [mul_inv_cancel_left]
  show f (Φ (fun i => g i * h i)) = f (Φ g)
  rw [key]
  exact e _ (hWV h hh)

/-- 核 —— 添字の付け替えは `IsNearOneEquivariant`。 -/
theorem isNearOneEquivariant_reindex (σ : ι → κ) :
    IsNearOneEquivariant (fun g : κ → G => fun i => g (σ i)) := by
  intro g V hV
  refine ⟨V, hV, fun h hh i => ?_⟩
  simpa using hh (σ i)

end AbstractCoreComp

section AbstractCoreContract

open Topology

variable {G : Type*} [Group G] [TopologicalSpace G] [IsTopologicalGroup G]

/-- ★★★核 —— **隣り合う 2 座標の積をとる写像は `IsNearOneEquivariant`**。

★これが `d` の第 2 項(`Fin.contractNth`)を扱う唯一の非自明な核である。
★使うのは位相群での共役の連続性(`ContinuousAt`)と `exists_nhds_one_split` だけ。 -/
theorem isNearOneEquivariant_contractNth {n : ℕ} (j : Fin (n + 1)) :
    IsNearOneEquivariant (fun g : Fin (n + 1) → G => Fin.contractNth j (· * ·) g) := by
  intro g V hV
  obtain ⟨V₁, hV₁, hsplit⟩ := exists_nhds_one_split hV
  refine ⟨(V ∩ V₁) ∩ ⋂ k : Fin n, (fun y => (g k.succ)⁻¹ * y * g k.succ) ⁻¹' V₁, ?_, ?_⟩
  · refine Filter.inter_mem (Filter.inter_mem hV hV₁) (Filter.iInter_mem.2 fun k => ?_)
    have hcont : ContinuousAt (fun y : G => (g k.succ)⁻¹ * y * g k.succ) 1 :=
      ((continuous_const.mul continuous_id).mul continuous_const).continuousAt
    exact hcont.preimage_mem_nhds (by simpa using hV₁)
  · intro h hh i
    show (Fin.contractNth j (· * ·) g i)⁻¹ *
      Fin.contractNth j (· * ·) (fun i => g i * h i) i ∈ V
    rcases lt_trichotomy (i : ℕ) (j : ℕ) with hlt | heq | hgt
    · rw [Fin.contractNth_apply_of_lt j _ g i hlt, Fin.contractNth_apply_of_lt j _ _ i hlt,
        inv_mul_cancel_left]
      exact (hh i.castSucc).1.1
    · rw [Fin.contractNth_apply_of_eq j _ g i heq, Fin.contractNth_apply_of_eq j _ _ i heq]
      have key : (g i.castSucc * g i.succ)⁻¹ * (g i.castSucc * h i.castSucc *
          (g i.succ * h i.succ))
          = ((g i.succ)⁻¹ * h i.castSucc * g i.succ) * h i.succ := by group
      rw [key]
      exact hsplit _ (Set.mem_iInter.1 (hh i.castSucc).2 i) _ (hh i.succ).1.2
    · rw [Fin.contractNth_apply_of_gt j _ g i hgt, Fin.contractNth_apply_of_gt j _ _ i hgt,
        inv_mul_cancel_left]
      exact (hh i.succ).1.1

end AbstractCoreContract

/-! ## 3. ★★★抽象核 III —— 作用の項(語彙ゼロ)

★`[Group G] [AddCommGroup M]` と `ρ : G → M →+ M` だけ(#209)。 -/

section AbstractCoreAction

open Topology

variable {G : Type*} [Group G] [TopologicalSpace G] {M : Type*} [AddCommGroup M]

/-- 核 —— 作用 `ρ` の**なめらかさ**: どの元も `1` のある近傍で固定される。
★`M` が離散のときの「作用の連続性」にほかならない(`smoothAction_iff_continuousAt`)。 -/
def SmoothAction (ρ : G → M →+ M) : Prop :=
  ∀ a : M, ∃ V ∈ 𝓝 (1 : G), ∀ u ∈ V, ρ u a = a

/-- ★★核 —— `d` の「作用の項」`g ↦ ρ (g c) (f (g ∘ σ))` はなめらかである。

★★点 `g` を固定すると `f (g ∘ σ)` は**ただ 1 つの元**なので、`SmoothAction` を
その元に使えば済む。★これが「一様性を要求しなくてよい」理由である。 -/
theorem LocConstNearOne.rhoTerm {ι κ : Type*} {ρ : G → M →+ M}
    (hmul : ∀ (x y : G) (a : M), ρ (x * y) a = ρ x (ρ y a)) (hρ : SmoothAction ρ)
    (c : κ) (σ : ι → κ) {f : (ι → G) → M} (hf : LocConstNearOne f) :
    LocConstNearOne (fun g : κ → G => ρ (g c) (f (fun i => g (σ i)))) := by
  intro g
  obtain ⟨V₁, hV₁, e₁⟩ := hf (fun i => g (σ i))
  obtain ⟨V₂, hV₂, e₂⟩ := hρ (f (fun i => g (σ i)))
  refine ⟨V₁ ∩ V₂, Filter.inter_mem hV₁ hV₂, fun h hh => ?_⟩
  have E := e₁ (fun i => h (σ i)) (fun i => (hh (σ i)).1)
  show ρ (g c * h c) (f (fun i => g (σ i) * h (σ i))) = ρ (g c) (f (fun i => g (σ i)))
  rw [E, hmul, e₂ (h c) (hh c).2]

end AbstractCoreAction

/-! ## 4. ★★★抽象核 IV —— 次数 1(語彙ゼロ)

★★**1-コサイクルについては「なめらか」と「`1` の近傍で消える」が同値**である。
★これが `CohomologyColimit.lean` の `VanishesNearOne` と本ファイルを繋ぐ。 -/

section AbstractCoreDegreeOne

open Topology

variable {G : Type*} [Group G] [TopologicalSpace G] {M : Type*} [AddCommGroup M]

omit [AddCommGroup M] in
/-- 核 —— 次数 1 のコチェインのなめらかさは「各点で右移動について局所定数」。 -/
theorem locConstNearOne_one_iff (f : G → M) :
    LocConstNearOne (fun g : Fin 1 → G => f (g 0))
      ↔ ∀ x : G, ∃ V ∈ 𝓝 (1 : G), ∀ v ∈ V, f (x * v) = f x := by
  constructor
  · intro hf x
    obtain ⟨V, hV, e⟩ := hf (fun _ => x)
    exact ⟨V, hV, fun v hv => e (fun _ => v) (fun _ => hv)⟩
  · intro hf g
    obtain ⟨V, hV, e⟩ := hf (g 0)
    exact ⟨V, hV, fun h hh => e (h 0) (hh 0)⟩

/-- ★★★核 —— **1-コサイクルについては、なめらかさは「`1` の近傍で消える」ことと同値**。

★仮定はコサイクル条件と `ρ 1 = id` だけ。★副有限性も離散性も使わない。 -/
theorem locConstNearOne_one_iff_vanishes {ρ : G → M →+ M} {f : G → M}
    (hρ1 : ∀ a, ρ 1 a = a) (hf : ∀ x y : G, f (x * y) = ρ x (f y) + f x) :
    LocConstNearOne (fun g : Fin 1 → G => f (g 0))
      ↔ ∃ V ∈ 𝓝 (1 : G), ∀ v ∈ V, f v = 0 := by
  rw [locConstNearOne_one_iff]
  have h11 := hf 1 1
  rw [one_mul, hρ1] at h11
  have hf1 : f 1 = 0 := (add_left_cancel (a := f 1) (b := 0) (c := f 1) (by
    rw [add_zero]; exact h11)).symm
  constructor
  · intro h
    obtain ⟨V, hV, hVf⟩ := h 1
    refine ⟨V, hV, fun v hv => ?_⟩
    have e := hVf v hv
    rw [one_mul] at e
    rw [e, hf1]
  · rintro ⟨V, hV, hVf⟩ x
    exact ⟨V, hV, fun v hv => by rw [hf x v, hVf v hv, map_zero, zero_add]⟩

/-- 核 —— 部分群の相対位相での「`1` の近傍で消える」の言い換え。
★`VanishesNearOne`(近傍は `G` の中)と、部分群 `↥S` の中の近傍を突き合わせる。 -/
theorem vanishesNearOne_iff_nhds_subtype (S : Subgroup G) (z : ↥S → M) :
    VanishesNearOne S z ↔ ∃ W ∈ 𝓝 (1 : ↥S), ∀ t ∈ W, z t = 0 := by
  constructor
  · rintro ⟨V, hV, hz⟩
    refine ⟨Subtype.val ⁻¹' V, ?_, fun t ht => hz t ht⟩
    rw [mem_nhds_subtype]
    exact ⟨V, by simpa using hV, le_rfl⟩
  · rintro ⟨W, hW, hz⟩
    rw [mem_nhds_subtype] at hW
    obtain ⟨V, hV, hVW⟩ := hW
    exact ⟨V, by simpa using hV, fun t ht => hz t (hVW ht)⟩

end AbstractCoreDegreeOne

/-! ## 5. ★★★抽象核 V —— 離散性(★`DiscreteTopology` を使う唯一の場所)

★★ここで `LocConstNearOne` が**「連続」そのもの**であることを両向き証明する。 -/

section DiscreteBridge

open Topology

variable {G : Type*} [Group G] [TopologicalSpace G] {M : Type*} [AddCommGroup M]
  {ι : Type*}

omit [AddCommGroup M] in
/-- ★★核 —— なめらかな関数は**連続**である(`M` は離散)。 -/
theorem continuous_of_locConstNearOne [Finite ι] [TopologicalSpace M] [DiscreteTopology M]
    [ContinuousMul G] {f : (ι → G) → M} (hf : LocConstNearOne f) : Continuous f := by
  rw [continuous_iff_continuousAt]
  intro g
  obtain ⟨V, hV, e⟩ := hf g
  have hU : (⋂ i, (fun x : ι → G => (g i)⁻¹ * x i) ⁻¹' V) ∈ 𝓝 g := by
    refine Filter.iInter_mem.2 fun i => ?_
    have hcont : ContinuousAt (fun x : ι → G => (g i)⁻¹ * x i) g :=
      (continuous_const.mul (continuous_apply i)).continuousAt
    exact hcont.preimage_mem_nhds (by simpa using hV)
  rw [ContinuousAt, nhds_discrete M, Filter.tendsto_pure]
  filter_upwards [hU] with x hx
  have hmem : ∀ i, (g i)⁻¹ * x i ∈ V := fun i => Set.mem_iInter.1 hx i
  have hx' : (fun i => g i * ((g i)⁻¹ * x i)) = x := by funext i; rw [mul_inv_cancel_left]
  rw [← hx']
  exact e _ hmem

omit [AddCommGroup M] in
/-- ★★核 —— 逆に、連続な関数はなめらかである(`M` は離散、添字は有限)。 -/
theorem locConstNearOne_of_continuous [Finite ι] [TopologicalSpace M] [DiscreteTopology M]
    [ContinuousMul G] {f : (ι → G) → M} (hf : Continuous f) : LocConstNearOne f := by
  intro g
  have hs : f ⁻¹' {f g} ∈ 𝓝 g :=
    hf.continuousAt.preimage_mem_nhds (by rw [nhds_discrete M]; exact rfl)
  rw [nhds_pi, Filter.mem_pi] at hs
  obtain ⟨I, -, t, ht, hsub⟩ := hs
  refine ⟨⋂ i, (fun y : G => g i * y) ⁻¹' (t i), Filter.iInter_mem.2 fun i => ?_, fun h hh => ?_⟩
  · have hcont : ContinuousAt (fun y : G => g i * y) 1 :=
      (continuous_const.mul continuous_id).continuousAt
    exact hcont.preimage_mem_nhds (by simpa using ht i)
  · have hpi : (fun i => g i * h i) ∈ Set.pi I t := fun i _ => Set.mem_iInter.1 (hh i) i
    simpa using hsub hpi

omit [AddCommGroup M] in
/-- ★★★核 —— **`LocConstNearOne` は「連続」そのもの**(`M` 離散・添字有限)。 -/
theorem locConstNearOne_iff_continuous [Finite ι] [TopologicalSpace M] [DiscreteTopology M]
    [ContinuousMul G] {f : (ι → G) → M} : LocConstNearOne f ↔ Continuous f :=
  ⟨continuous_of_locConstNearOne, locConstNearOne_of_continuous⟩

/-- ★核 —— `SmoothAction` は「作用が `1` で連続」と同値(`M` 離散)。 -/
theorem smoothAction_iff_continuousAt [TopologicalSpace M] [DiscreteTopology M] {ρ : G → M →+ M}
    (hρ1 : ∀ a, ρ 1 a = a) : SmoothAction ρ ↔ ∀ a : M, ContinuousAt (fun g : G => ρ g a) 1 := by
  constructor
  · rintro h a
    obtain ⟨V, hV, hVa⟩ := h a
    rw [ContinuousAt, hρ1, nhds_discrete M, Filter.tendsto_pure]
    filter_upwards [hV] with u hu using hVa u hu
  · intro h a
    have hs : (fun g : G => ρ g a) ⁻¹' {a} ∈ 𝓝 (1 : G) := by
      refine (h a).preimage_mem_nhds ?_
      rw [hρ1, nhds_discrete M]
      exact rfl
    exact ⟨_, hs, fun u hu => hu⟩

end DiscreteBridge

/-! ## 6. ★★★★具体層 I —— 連続コチェイン `C^n_cont(G,A)`

★★ここから `Rep k G` に代入する。★宇宙は `Rep.{0} k G` に固定する(逸脱 3、#224)。 -/

section ContinuousCochains

open Topology

variable {k G : Type} [CommRing k] [Group G] [TopologicalSpace G]

/-- ★★★★**連続コチェイン** `C^n_cont(G,A) ⊆ C^n(G,A)`。

★`LocConstNearOne` で定義するので、係数 `A` の位相は要求しない。
★`↥A` が離散なら `Continuous` と同値(`mem_contCochains_iff_continuous`)。 -/
def contCochains (A : Rep.{0} k G) (n : ℕ) : Submodule k ((Fin n → G) → ↥A) :=
  smoothFuns k (Fin n)

theorem mem_contCochains {A : Rep.{0} k G} {n : ℕ} {f : (Fin n → G) → ↥A} :
    f ∈ contCochains A n ↔ LocConstNearOne f := Iff.rfl

omit [TopologicalSpace G] in
/-- コバウンダリ作用素の値(★`rfl`)。 -/
theorem inhomogeneousCochains_d_apply (A : Rep.{0} k G) (n : ℕ) (f : (Fin n → G) → ↥A)
    (g : Fin (n + 1) → G) :
    ModuleCat.Hom.hom (inhomogeneousCochains.d A n) f g
      = A.ρ (g 0) (f fun i => g i.succ)
        + ∑ j : Fin (n + 1), (-1 : k) ^ ((j : ℕ) + 1) • f (Fin.contractNth j (· * ·) g) := rfl

omit [TopologicalSpace G] in
/-- コバウンダリ作用素の 2 つの項への分解。 -/
theorem inhomogeneousCochains_d_eq (A : Rep.{0} k G) (n : ℕ) (f : (Fin n → G) → ↥A) :
    ModuleCat.Hom.hom (inhomogeneousCochains.d A n) f
      = (fun g : Fin (n + 1) → G => repAddHom A (g 0) (f fun i => g i.succ))
        + ∑ j : Fin (n + 1), ((-1 : k) ^ ((j : ℕ) + 1) •
            (fun g : Fin (n + 1) → G => f (Fin.contractNth j (· * ·) g))) := by
  funext g
  simp only [Pi.add_apply, Finset.sum_apply, Pi.smul_apply]
  rfl

/-- ★★★★★**コバウンダリ作用素は連続コチェインを連続コチェインに写す**。

★★これが本ファイルの中核である。★中身は 3 本の抽象核
(`LocConstNearOne.rhoTerm` / `LocConstNearOne.comp` / `isNearOneEquivariant_contractNth`)
に入っており、ここでは代入するだけ。 -/
theorem d_mem_contCochains [IsTopologicalGroup G] (A : Rep.{0} k G)
    (hA : SmoothAction (repAddHom A)) (n : ℕ) {f : (Fin n → G) → ↥A}
    (hf : f ∈ contCochains A n) :
    ModuleCat.Hom.hom (inhomogeneousCochains.d A n) f ∈ contCochains A (n + 1) := by
  rw [inhomogeneousCochains_d_eq]
  refine Submodule.add_mem _ ?_ (Submodule.sum_mem _ fun j _ => Submodule.smul_mem _ _ ?_)
  · exact LocConstNearOne.rhoTerm (repAddHom_mul A) hA 0 Fin.succ hf
  · exact hf.comp (isNearOneEquivariant_contractNth j)

/-- ★退化の自己検査 —— `C^0_cont(G,A) = A`(次数 0 に条件は無い)。
★★「なめらかベクトル」に痩せていないことの確認である。 -/
theorem contCochains_zero (A : Rep.{0} k G) : contCochains A 0 = ⊤ := by
  refine eq_top_iff.2 fun f _ g => ⟨Set.univ, Filter.univ_mem, fun h _ => ?_⟩
  exact congrArg f (Subsingleton.elim _ _)

/-- ★退化の自己検査 —— `G` が離散なら連続コチェインはすべてのコチェイン。 -/
theorem contCochains_of_discreteTopology (A : Rep.{0} k G) [DiscreteTopology G] (n : ℕ) :
    contCochains A n = ⊤ := by
  refine eq_top_iff.2 fun f _ g => ⟨{1}, ?_, fun h hh => ?_⟩
  · exact (isOpen_discrete ({1} : Set G)).mem_nhds rfl
  · have hg : (fun i => g i * h i) = g := by
      funext i; rw [hh i, mul_one]
    rw [hg]

/-- ★★★★**`C^n_cont(G,A)` は本当に「連続コチェイン」である**
(★`↥A` の離散性を使う具体層)。 -/
theorem mem_contCochains_iff_continuous (A : Rep.{0} k G) [TopologicalSpace ↥A]
    [DiscreteTopology ↥A] [ContinuousMul G] (n : ℕ) (f : (Fin n → G) → ↥A) :
    f ∈ contCochains A n ↔ Continuous f :=
  locConstNearOne_iff_continuous

/-- ★★`SmoothAction (repAddHom A)` は「`G` の `A` への作用が連続」と同値(`↥A` 離散)。 -/
theorem smoothAction_repAddHom_iff (A : Rep.{0} k G) [TopologicalSpace ↥A]
    [DiscreteTopology ↥A] :
    SmoothAction (repAddHom A) ↔ ∀ a : ↥A, ContinuousAt (fun g : G => A.ρ g a) 1 :=
  smoothAction_iff_continuousAt (repAddHom_one A)

end ContinuousCochains

/-! ## 7. ★★★★★具体層 II —— `Z^n_cont` / `B^n_cont` / `H^n_cont` -/

section ContinuousCohomology

open Topology

variable {k G : Type} [CommRing k] [Group G] [TopologicalSpace G]

/-- ★★★**連続コサイクル** `Z^n_cont(G,A)`。 -/
def contCocycles (A : Rep.{0} k G) (n : ℕ) : Submodule k ((Fin n → G) → ↥A) :=
  contCochains A n ⊓ LinearMap.ker (ModuleCat.Hom.hom (inhomogeneousCochains.d A n))

theorem mem_contCocycles {A : Rep.{0} k G} {n : ℕ} {f : (Fin n → G) → ↥A} :
    f ∈ contCocycles A n ↔
      LocConstNearOne f ∧ ModuleCat.Hom.hom (inhomogeneousCochains.d A n) f = 0 := Iff.rfl

/-- ★★★**連続コバウンダリ** —— `Z^n_cont(G,A)` の部分加群として。

★★`d` が連続コチェインを保つ仮定(`SmoothAction`)を**定義には課さない**。
仮定つきの正しい記述は `range_contD` である。 -/
def contCoboundaries (A : Rep.{0} k G) : (n : ℕ) → Submodule k ↥(contCocycles A n)
  | 0 => ⊥
  | n + 1 => Submodule.comap (contCocycles A (n + 1)).subtype
      (Submodule.map (ModuleCat.Hom.hom (inhomogeneousCochains.d A n)) (contCochains A n))

@[simp] theorem contCoboundaries_zero (A : Rep.{0} k G) : contCoboundaries A 0 = ⊥ := rfl

theorem contCoboundaries_succ (A : Rep.{0} k G) (n : ℕ) :
    contCoboundaries A (n + 1) = Submodule.comap (contCocycles A (n + 1)).subtype
      (Submodule.map (ModuleCat.Hom.hom (inhomogeneousCochains.d A n)) (contCochains A n)) := rfl

/-- ★★★★★**連続コホモロジー** `H^n_cont(G,A) = Z^n_cont ⧸ B^n_cont`。 -/
def Hcont (A : Rep.{0} k G) (n : ℕ) : Type :=
  ↥(contCocycles A n) ⧸ contCoboundaries A n

instance (A : Rep.{0} k G) (n : ℕ) : AddCommGroup (Hcont A n) :=
  inferInstanceAs (AddCommGroup (_ ⧸ contCoboundaries A n))

instance (A : Rep.{0} k G) (n : ℕ) : Module k (Hcont A n) :=
  inferInstanceAs (Module k (_ ⧸ contCoboundaries A n))

/-- 連続コサイクルの類。 -/
def Hcontπ (A : Rep.{0} k G) (n : ℕ) : ↥(contCocycles A n) →ₗ[k] Hcont A n :=
  (contCoboundaries A n).mkQ

theorem Hcontπ_surjective (A : Rep.{0} k G) (n : ℕ) : Function.Surjective (Hcontπ A n) :=
  Submodule.mkQ_surjective _

omit [TopologicalSpace G] in
/-- ★`d ∘ d = 0`(mathlib の非斉次コチェイン複体から)。 -/
theorem d_comp_d_apply (A : Rep.{0} k G) (n : ℕ) (f : (Fin n → G) → ↥A) :
    ModuleCat.Hom.hom (inhomogeneousCochains.d A (n + 1))
      (ModuleCat.Hom.hom (inhomogeneousCochains.d A n) f) = 0 := by
  have h := groupCohomology.inhomogeneousCochains.d_comp_d (A := A) (n := n)
  have h2 := congrArg (fun (φ : ModuleCat.of k ((Fin n → G) → (A : Type)) ⟶ _) =>
    ModuleCat.Hom.hom φ f) h
  simpa using h2

/-- ★★連続コチェインのコバウンダリは**連続コサイクル**である。 -/
theorem d_mem_contCocycles [IsTopologicalGroup G] (A : Rep.{0} k G)
    (hA : SmoothAction (repAddHom A)) (n : ℕ) {f : (Fin n → G) → ↥A}
    (hf : f ∈ contCochains A n) :
    ModuleCat.Hom.hom (inhomogeneousCochains.d A n) f ∈ contCocycles A (n + 1) :=
  ⟨d_mem_contCochains A hA n hf, LinearMap.mem_ker.2 (d_comp_d_apply A n f)⟩

theorem mem_contCoboundaries_succ_iff (A : Rep.{0} k G) (n : ℕ) (z : ↥(contCocycles A (n + 1))) :
    z ∈ contCoboundaries A (n + 1) ↔ ∃ f ∈ contCochains A n,
      ModuleCat.Hom.hom (inhomogeneousCochains.d A n) f = (z : (Fin (n + 1) → G) → ↥A) := by
  rw [contCoboundaries_succ, Submodule.mem_comap]
  exact Submodule.mem_map

/-- ★連続コチェインの上のコバウンダリ写像 `C^n_cont → Z^{n+1}_cont`。 -/
noncomputable def contD [IsTopologicalGroup G] (A : Rep.{0} k G)
    (hA : SmoothAction (repAddHom A)) (n : ℕ) :
    ↥(contCochains A n) →ₗ[k] ↥(contCocycles A (n + 1)) :=
  LinearMap.restrict (ModuleCat.Hom.hom (inhomogeneousCochains.d A n))
    (fun _ hx => d_mem_contCocycles A hA n hx)

/-- ★★`B^{n+1}_cont = im(d : C^n_cont → Z^{n+1}_cont)`。 -/
theorem range_contD [IsTopologicalGroup G] (A : Rep.{0} k G)
    (hA : SmoothAction (repAddHom A)) (n : ℕ) :
    LinearMap.range (contD A hA n) = contCoboundaries A (n + 1) := by
  ext z
  rw [LinearMap.mem_range, mem_contCoboundaries_succ_iff]
  constructor
  · rintro ⟨⟨f, hf⟩, rfl⟩
    exact ⟨f, hf, rfl⟩
  · rintro ⟨f, hf, hfz⟩
    exact ⟨⟨f, hf⟩, Subtype.ext hfz⟩

/-- ★★★退化の自己検査 —— **`Z^0_cont(G,A) = A^G`**。 -/
theorem mem_contCocycles_zero_iff (A : Rep.{0} k G) (f : (Fin 0 → G) → ↥A) :
    f ∈ contCocycles A 0 ↔ ∀ g : G, A.ρ g (f default) = f default := by
  have hsub : ∀ x : Fin 0 → G, f x = f default := fun x => congrArg f (Subsingleton.elim x default)
  constructor
  · rintro ⟨-, hker⟩ g
    have h : ModuleCat.Hom.hom (inhomogeneousCochains.d A 0) f (fun _ => g)
        = (0 : (Fin (0 + 1) → G) → ↥A) (fun _ => g) := congrFun hker _
    rw [inhomogeneousCochains_d_apply] at h
    simp only [hsub, Pi.zero_apply] at h
    simpa [add_neg_eq_zero] using h
  · intro h
    refine ⟨contCochains_zero A ▸ Submodule.mem_top, LinearMap.mem_ker.2 ?_⟩
    show ModuleCat.Hom.hom (inhomogeneousCochains.d A 0) f = (0 : (Fin (0 + 1) → G) → ↥A)
    funext g
    rw [inhomogeneousCochains_d_apply]
    simp only [hsub, Pi.zero_apply]
    simp [h (g 0)]

/-- ★★★退化の自己検査 —— **`H^0_cont(G,A) = Z^0_cont(G,A) = A^G`**。 -/
noncomputable def HcontZeroEquiv (A : Rep.{0} k G) : Hcont A 0 ≃ₗ[k] ↥(contCocycles A 0) :=
  Submodule.quotEquivOfEqBot _ (contCoboundaries_zero A)

end ContinuousCohomology

/-! ## 8. ★★★★★具体層 III —— 次数 1 に潰すと `smoothCocycles₁` に一致する -/

section DegreeOne

open Topology

variable {k H : Type} [CommRing k] [Group H] [TopologicalSpace H]

omit [TopologicalSpace H] in
/-- 次数 1 でのコバウンダリ作用素の値。 -/
theorem d_one_apply (B : Rep.{0} k H) (F : (Fin 1 → H) → ↥B) (g : Fin 2 → H) :
    ModuleCat.Hom.hom (inhomogeneousCochains.d B 1) F g
      = B.ρ (g 0) (F fun _ => g 1) - F (fun _ => g 0 * g 1) + F (fun _ => g 0) := by
  have e0 : (fun i : Fin 1 => g i.succ) = fun _ => g 1 := by
    funext i; rw [Subsingleton.elim i 0]; rfl
  have e1 : Fin.contractNth (0 : Fin 2) (· * ·) g = fun _ => g 0 * g 1 := by
    funext j
    rw [Fin.contractNth_apply_of_eq _ _ _ j (by rw [Subsingleton.elim j 0]; rfl)]
    rw [Subsingleton.elim j 0]; rfl
  have e2 : Fin.contractNth (1 : Fin 2) (· * ·) g = fun _ => g 0 := by
    funext j
    rw [Fin.contractNth_apply_of_lt _ _ _ j (by rw [Subsingleton.elim j 0]; norm_num)]
    rw [Subsingleton.elim j 0]; rfl
  rw [inhomogeneousCochains_d_apply, Fin.sum_univ_two, e0, e1, e2]
  norm_num
  abel

omit [TopologicalSpace H] in
/-- ★★次数 1 のコサイクル条件は mathlib の `cocycles₁` と一致する。 -/
theorem mem_ker_d_one_iff (B : Rep.{0} k H) (F : (Fin 1 → H) → ↥B) :
    F ∈ LinearMap.ker (ModuleCat.Hom.hom (inhomogeneousCochains.d B 1))
      ↔ (fun h : H => F fun _ => h) ∈ cocycles₁ B := by
  rw [mem_cocycles₁_iff]
  constructor
  · intro hker x y
    have h : ModuleCat.Hom.hom (inhomogeneousCochains.d B 1) F ![x, y]
        = (0 : (Fin 2 → H) → ↥B) ![x, y] := congrFun (LinearMap.mem_ker.1 hker) _
    rw [d_one_apply] at h
    simp only [Matrix.cons_val_zero, Matrix.cons_val_one, Pi.zero_apply] at h
    rw [sub_add_eq_add_sub, sub_eq_zero] at h
    exact h.symm
  · intro hc
    refine LinearMap.mem_ker.2 ?_
    show ModuleCat.Hom.hom (inhomogeneousCochains.d B 1) F = (0 : (Fin 2 → H) → ↥B)
    funext g
    have hxy : F (fun _ => g 0 * g 1) = B.ρ (g 0) (F fun _ => g 1) + F (fun _ => g 0) :=
      hc (g 0) (g 1)
    rw [d_one_apply, hxy]
    simp

end DegreeOne

section DegreeOneCompare

open Topology

variable {k G : Type} [CommRing k] [Group G] [TopologicalSpace G]

/-- ★★★★★**退化の自己検査 —— 次数 1 に潰すと `smoothCocycles₁` に一致する**。

左辺は本ファイルの連続コチェイン(位相群 `↥S` の上、次数 1)、
右辺は `CohomologyColimit.lean` の「`1` の近傍で消える 1-コサイクル」。
★★仮定は `z` がコサイクルであることだけ。副有限性も離散性も要らない。 -/
theorem mem_contCochains_one_iff_mem_smoothCocycles₁ (A : Rep.{0} k G) (S : Subgroup G)
    (z : cocycles₁ (Rep.res S.subtype A)) :
    (fun g : Fin 1 → ↥S => (z : ↥S → ↥A) (g 0)) ∈ contCochains (Rep.res S.subtype A) 1
      ↔ z ∈ smoothCocycles₁ A S := by
  have h := locConstNearOne_one_iff_vanishes (ρ := repAddHom (Rep.res S.subtype A))
    (f := fun t : ↥S => (z : ↥S → ↥A) t) (repAddHom_one _)
    (fun x y => (mem_cocycles₁_iff (A := Rep.res S.subtype A) _).1 z.2 x y)
  exact h.trans (vanishesNearOne_iff_nhds_subtype S ((z : ↥S → ↥A))).symm

end DegreeOneCompare

/-! ## 9. ★★★★具体層 IV —— inflation(colimit 側との比較射)

★★**`S` が `1` の近傍なら inflation は連続コチェインを与える**。
★これが `colim_i H^n(G ⧸ S i, A^{S i}) → H^n_cont(G,A)` が存在する理由である。 -/

section Inflation

open Topology

variable {k G : Type} [CommRing k] [Group G] [TopologicalSpace G]

/-- **inflation のコチェイン水準** `C^n(G ⧸ S, A^S) → C^n(G,A)`。 -/
noncomputable def inflCochain (A : Rep.{0} k G) (S : Subgroup G) [S.Normal] (n : ℕ) :
    ((Fin n → G ⧸ S) → ↥(A.quotientToInvariants S)) →ₗ[k] ((Fin n → G) → ↥A) :=
  ModuleCat.Hom.hom ((cochainsMap (A := A.quotientToInvariants S) (B := A) (QuotientGroup.mk' S)
    (Rep.ofHom (A.ρ.quotientToInvariants_lift S))).f n)

omit [TopologicalSpace G] in
theorem inflCochain_apply (A : Rep.{0} k G) (S : Subgroup G) [S.Normal] (n : ℕ)
    (q : (Fin n → G ⧸ S) → ↥(A.quotientToInvariants S)) (g : Fin n → G) :
    inflCochain A S n q g
      = (A.ρ.quotientToInvariants_lift S) (q fun i => QuotientGroup.mk' S (g i)) := rfl

omit [TopologicalSpace G] in
/-- ★★inflation はコバウンダリ作用素と可換(★mathlib の `cochainsMap` の可換性)。 -/
theorem inflCochain_comm_d (A : Rep.{0} k G) (S : Subgroup G) [S.Normal] (n : ℕ)
    (q : (Fin n → G ⧸ S) → ↥(A.quotientToInvariants S)) :
    ModuleCat.Hom.hom (inhomogeneousCochains.d A n) (inflCochain A S n q)
      = inflCochain A S (n + 1)
          (ModuleCat.Hom.hom (inhomogeneousCochains.d (A.quotientToInvariants S) n) q) := by
  have h := (cochainsMap (A := A.quotientToInvariants S) (B := A) (QuotientGroup.mk' S)
    (Rep.ofHom (A.ρ.quotientToInvariants_lift S))).comm n (n + 1)
  rw [inhomogeneousCochains.d_def, inhomogeneousCochains.d_def] at h
  have h2 := congrArg (fun φ => ModuleCat.Hom.hom φ q) h
  simp only [ModuleCat.hom_comp, LinearMap.comp_apply] at h2
  exact h2

/-- ★★★★**`S` が `1` の近傍なら、inflation したコチェインは連続である**。

★★これが (c)-3 の「比較関手」の中身である。★開性も副有限性も要らず、
`(S : Set G) ∈ 𝓝 1` だけで足りる。 -/
theorem inflCochain_mem_contCochains (A : Rep.{0} k G) (S : Subgroup G) [S.Normal]
    (hS : (S : Set G) ∈ 𝓝 (1 : G)) (n : ℕ)
    (q : (Fin n → G ⧸ S) → ↥(A.quotientToInvariants S)) :
    inflCochain A S n q ∈ contCochains A n := by
  intro g
  refine ⟨(S : Set G), hS, fun h hh => ?_⟩
  rw [inflCochain_apply, inflCochain_apply]
  congr 2
  funext i
  show QuotientGroup.mk' S (g i * h i) = QuotientGroup.mk' S (g i)
  rw [map_mul, QuotientGroup.mk'_apply, QuotientGroup.mk'_apply,
    (QuotientGroup.eq_one_iff (h i)).2 (hh i), mul_one]

/-- ★★★連続コサイクルへの inflation。 -/
theorem inflCochain_mem_contCocycles (A : Rep.{0} k G) (S : Subgroup G) [S.Normal]
    (hS : (S : Set G) ∈ 𝓝 (1 : G)) (n : ℕ)
    (q : (Fin n → G ⧸ S) → ↥(A.quotientToInvariants S))
    (hq : ModuleCat.Hom.hom (inhomogeneousCochains.d (A.quotientToInvariants S) n) q = 0) :
    inflCochain A S n q ∈ contCocycles A n :=
  ⟨inflCochain_mem_contCochains A S hS n q,
    LinearMap.mem_ker.2 (by rw [inflCochain_comm_d, hq, map_zero])⟩

/-- ★★inflation はコバウンダリを**連続**コバウンダリに写す。 -/
theorem inflCochain_d_mem_contCoboundaries (A : Rep.{0} k G) (S : Subgroup G) [S.Normal]
    (hS : (S : Set G) ∈ 𝓝 (1 : G)) (n : ℕ)
    (p : (Fin n → G ⧸ S) → ↥(A.quotientToInvariants S)) (z : ↥(contCocycles A (n + 1)))
    (hzp : (z : (Fin (n + 1) → G) → ↥A)
      = inflCochain A S (n + 1)
          (ModuleCat.Hom.hom (inhomogeneousCochains.d (A.quotientToInvariants S) n) p)) :
    z ∈ contCoboundaries A (n + 1) := by
  rw [mem_contCoboundaries_succ_iff]
  exact ⟨inflCochain A S n p, inflCochain_mem_contCochains A S hS n p,
    by rw [inflCochain_comm_d, hzp]⟩

/-- 塔の 1 段の inflation(コチェイン水準) `C^n(G ⧸ S, A^S) → C^n(G ⧸ T, A^T)`。 -/
noncomputable def inflStepCochain (A : Rep.{0} k G) (S T : Subgroup G) [S.Normal] [T.Normal]
    (hle : T ≤ S) (n : ℕ) :
    ((Fin n → G ⧸ S) → ↥(A.quotientToInvariants S)) →ₗ[k]
      ((Fin n → G ⧸ T) → ↥(A.quotientToInvariants T)) :=
  ModuleCat.Hom.hom ((cochainsMap (A := A.quotientToInvariants S) (B := A.quotientToInvariants T)
    (quotStep S T hle) (invStep A S T hle)).f n)

omit [TopologicalSpace G] in
/-- ★★★**塔と両立する** —— `T ≤ S` のとき `G ⧸ T` を経由した inflation は
`G ⧸ S` からの inflation に一致する。★**証明は `rfl`**。 -/
theorem inflCochain_inflStepCochain (A : Rep.{0} k G) (S T : Subgroup G) [S.Normal] [T.Normal]
    (hle : T ≤ S) (n : ℕ) (q : (Fin n → G ⧸ S) → ↥(A.quotientToInvariants S)) :
    inflCochain A T n (inflStepCochain A S T hle n q) = inflCochain A S n q := rfl

omit [TopologicalSpace G] in
/-- 塔の 1 段の inflation もコバウンダリ作用素と可換。 -/
theorem inflStepCochain_comm_d (A : Rep.{0} k G) (S T : Subgroup G) [S.Normal] [T.Normal]
    (hle : T ≤ S) (n : ℕ) (q : (Fin n → G ⧸ S) → ↥(A.quotientToInvariants S)) :
    ModuleCat.Hom.hom (inhomogeneousCochains.d (A.quotientToInvariants T) n)
        (inflStepCochain A S T hle n q)
      = inflStepCochain A S T hle (n + 1)
          (ModuleCat.Hom.hom (inhomogeneousCochains.d (A.quotientToInvariants S) n) q) := by
  have h := (cochainsMap (A := A.quotientToInvariants S) (B := A.quotientToInvariants T)
    (quotStep S T hle) (invStep A S T hle)).comm n (n + 1)
  rw [inhomogeneousCochains.d_def, inhomogeneousCochains.d_def] at h
  have h2 := congrArg (fun φ => ModuleCat.Hom.hom φ q) h
  simp only [ModuleCat.hom_comp, LinearMap.comp_apply] at h2
  exact h2

/-- ★★★★★**比較射(コサイクルの水準)** `Z^n(G ⧸ S, A^S) → H^n_cont(G,A)`。 -/
noncomputable def inflHcont (A : Rep.{0} k G) (S : Subgroup G) [S.Normal]
    (hS : (S : Set G) ∈ 𝓝 (1 : G)) (n : ℕ) :
    ↥(LinearMap.ker (ModuleCat.Hom.hom (inhomogeneousCochains.d (A.quotientToInvariants S) n)))
      →ₗ[k] Hcont A n :=
  (Hcontπ A n).comp (LinearMap.restrict (inflCochain A S n)
    (fun q hq => inflCochain_mem_contCocycles A S hS n q (LinearMap.mem_ker.1 hq)))

/-- ★★★**比較射はコバウンダリを殺す** —— したがって `H^n(G ⧸ S, A^S)` の上で定まる。 -/
theorem inflHcont_d_eq_zero (A : Rep.{0} k G) (S : Subgroup G) [S.Normal]
    (hS : (S : Set G) ∈ 𝓝 (1 : G)) (n : ℕ)
    (p : (Fin n → G ⧸ S) → ↥(A.quotientToInvariants S)) :
    inflHcont A S hS (n + 1)
        ⟨ModuleCat.Hom.hom (inhomogeneousCochains.d (A.quotientToInvariants S) n) p,
          LinearMap.mem_ker.2 (d_comp_d_apply _ n p)⟩ = 0 := by
  refine (Submodule.Quotient.mk_eq_zero _).2 ?_
  exact inflCochain_d_mem_contCoboundaries A S hS n p _ rfl

/-- ★★★★**比較射は塔と両立する** —— 有向系(`CohomologyColimit.lean` の `h2Sys`)の
1 段を通しても値は変わらない。
★これが `colim_i H^n(G ⧸ S i, A^{S i}) → H^n_cont(G,A)` が定まる理由である。 -/
theorem inflHcont_inflStepCochain (A : Rep.{0} k G) (S T : Subgroup G) [S.Normal] [T.Normal]
    (hle : T ≤ S) (hS : (S : Set G) ∈ 𝓝 (1 : G)) (hT : (T : Set G) ∈ 𝓝 (1 : G)) (n : ℕ)
    (z : ↥(LinearMap.ker
      (ModuleCat.Hom.hom (inhomogeneousCochains.d (A.quotientToInvariants S) n)))) :
    inflHcont A T hT n ⟨inflStepCochain A S T hle n (z : _), by
        rw [LinearMap.mem_ker, inflStepCochain_comm_d, LinearMap.mem_ker.1 z.2, map_zero]⟩
      = inflHcont A S hS n z := by
  show Hcontπ A n _ = Hcontπ A n _
  exact congrArg (Hcontπ A n) (Subtype.ext (inflCochain_inflStepCochain A S T hle n (z : _)))

end Inflation

/-! ## 10. ★副有限の塔との接続(`CohomologyColimit.lean` の実例) -/

section ProfiniteTower

open Topology

variable {G : Type} [Group G] [TopologicalSpace G]

/-- 開部分群は `1` の近傍である。 -/
theorem coe_mem_nhds_one (S : Subgroup G) (hS : IsOpen (S : Set G)) :
    (S : Set G) ∈ 𝓝 (1 : G) := hS.mem_nhds S.one_mem

/-- ★`CohomologyColimit.lean` の副有限の塔の各段は `1` の近傍である。
★したがって §9 の inflation はその塔にそのまま適用できる。 -/
theorem profiniteTower_mem_nhds_one (H : (OpenNormalSubgroup G)ᵒᵈ) :
    ((profiniteTower G H : Subgroup G) : Set G) ∈ 𝓝 (1 : G) :=
  coe_mem_nhds_one _ (OrderDual.ofDual H : OpenNormalSubgroup G).isOpen'

end ProfiniteTower

end ABC3.Found.PGC
