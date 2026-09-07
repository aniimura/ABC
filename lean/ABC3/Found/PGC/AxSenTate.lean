import ABC3.Found.PGC.Section3RealParameters

/-!
# [pGC] Ax–Sen–Tate と Hodge-Tate 重み空間の有限次元性

`Found/PGC/Section3RealParameters.lean` は `hodgeTateWeightSpace` / `hodgeTateDim` /
`IsHodgeTate` を実物(自由パラメータ 0)で書いたが、その docstring が
自分で名指ししたとおり **`d_V(i)` の値については何も言えていなかった**:

> ★★`d_V(i)` の値については何も言えていない。 `d_{triv}(0) = 1` すら
> `ℂ_K^{Γ_K} = K`(Ax–Sen–Tate)を要し、それは mathlib にも木にも無い。

本ファイルはその穴を **1 つの仮説 `AxSenTate` に閉じ込め**、
「その仮説だけが足りない」ことを**証明**する。

## 在庫調査(自分で測った。コマンドを残す)

* mathlib: **無い**。
  `grep -i -E "AxSen|Ax_Sen|SenTate|axSen" .cache/mathlib-index.txt` → **0 件**。
  `grep -i "hodgeTate" .cache/mathlib-index.txt` → **0 件**。
  周辺(`RingTheory/Perfectoid/`)には `WittVector.fontaineTheta` /
  `BDeRham` があるが、`ℂ_K^{Γ_K} = K` に相当するものは無い。
* 木: **無い**。`grep -rn "AxSen" lean/ABC3/` は
  `ClosureCompletion.lean` の 2 つの **docstring** と
  `Section3RealParameters.lean` の「無い」という記述だけ。
* `Check/PGC/`(何が偽と分かっているかの台帳)にも Ax–Sen–Tate は無い
  ——本ファイルの主張は台帳の反例と衝突しない。

## 何が言えたか

### §1 抽象核(分岐・付値・Galois・p 進の語彙が **1 つも出てこない**)

`linearIndependent_of_fixedRange` ——
体 `A` に「型 `G` で添字された環準同型の族 `aut`」が作用し、`A`-加群 `M` に
加法準同型の族 `act` が**半線型**に作用しているとする。

* `hfix`: `aut` で固定される `A` の元は部分体 `F` から来る、
* `hw`: `w i` はどれも**共通の固有値** `u g ≠ 0` の固有ベクトル

ならば **`F`-一次独立 ⇒ `A`-一次独立**。

★`G` に群構造は要らない(`act` の乗法性も `aut` の単射性も使わない)。
★証明は「関係式の台を最小に取り、`σ` を掛けて引く」という古典的な議論を
`Finset.card` についての強帰納法で書いたもの。

`rank_le_finrank_of_fixedRange` / `finrank_le_of_fixedRange` /
`finite_of_fixedRange` はその系(`Module.Finite A M` のとき次元が押さえられる)。

★**非空虚性**: `linearIndependent_complex_of_real` ——
`F = ℝ`, `A = ℂ`, `aut = 複素共役`, `M = ι → ℂ` で抽象核の仮説は**全部真になる**
(結論は「実ベクトルが `ℝ` 上独立なら `ℂ` 上も独立」という古典的事実)。
★仮説が空虚でないことの証拠である。

### §2 具体層

* `AxSenTate K` —— `ℂ_K^{Γ_K} ⊆ K`。**逆向き `K ⊆ ℂ_K^{Γ_K}` は無条件に真**
  (`algebraMap_mem_fixedSubmodule`)なので、これは `ℂ_K^{Γ_K} = K` と同値。
* ★★`hodgeTateDim_trivial_zero_eq_one_iff` ——
  **`d_triv(0) = 1 ↔ AxSenTate K`**。
  ★★これが本ファイルの主結果である。「`d_V(i)` について何も言えない」原因が
  **ちょうど Ax–Sen–Tate ただ 1 つ**であって、それ以外に何も足りていないことを
  `↔` で示している(仮説が**弱すぎも強すぎもしない**)。
* ★`hodgeTateDim_le_finrank` / `finiteDimensional_hodgeTateWeightSpace` ——
  `AxSenTate K` を仮定すれば、有限次元 `V` について
  **`d_V(i) ≤ dim_{ℚ_p} V`** かつ重み空間は**有限次元**。
  ★`Module.finrank` が無限次元で `0` を返す問題はこれで解消する。
* ★`isHodgeTate_trivial` —— `AxSenTate K` の下で自明表現は Hodge-Tate。
  ★`IsHodgeTate` を**実際に満たす表現の最初の 1 つ**である
  (`Section3RealParameters.lean` は定義しただけだった)。
* ★`fixedSubmodule_ne_bot` —— `ℂ_K^{Γ_K} ≠ 0`(具体層の非空虚性)。

## ★原典より短い道(自分で見つけた点。★誇張しないように書く)

`d_V(i) ≤ dim V` そのものは古典的にも Ax–Sen–Tate だけから出る議論なので、
「原典より短い」と言えるのは**書き方**の方である:

1. ★**`G` に群構造が要らない。** 原典は `Γ_K` の Galois 群としての性質を通して
   議論を書くが、実際に使うのは「`act g` が加法的」「`aut g` が環準同型」の 2 つだけで、
   **合成も逆元も単位元も使わない**。`linearIndependent_of_fixedRange` の `G` は
   `Type*` のままである。
2. ★**`aut g` の単射性も使わない。** 台が縮むことは `aut g 0 = 0` だけで出る。
3. ★**Galois コホモロジーを 1 行も書いていない。** Hilbert 90
   (`H^1(Γ_K, ℂ_K^×) = 1`)も Tate の `H^i(Γ_K, ℂ_K(j))` の計算も使っていない。
   必要なのは「一次独立性の降下」だけで、それは **`Finset.card` の強帰納法**で閉じる。
4. ★**`Finsupp` を経由しない。** `linearIndependent_iff''`(`Finset` と `ι → A` の形)で
   書くと `Finsupp.mapRange` / `Finsupp.sum_mapRange_index` が 1 つも要らなくなる。

★逆に「短くならなかった」点も書いておく: `d_triv(0) = 1` は
**Ax–Sen–Tate と同値**であって(`hodgeTateDim_trivial_zero_eq_one_iff`)、
これを迂回する軽い道は**無い**。★本体の見立て(「重み空間の有限次元性の方が Ax–Sen–Tate より軽い」)は
**外れた** —— 有限次元性の証明も Ax–Sen–Tate を仮説として要求する
(`hodgeTateDim_le_finrank` の `hK`)。

## 逸脱の記録(CLAUDE.md「逸脱」)

1. **`AxSenTate` は仮説であって定理ではない。** 本ファイルは証明していない。
   証明には Ax の補題(「ほとんど不変な元はほとんど `K` に居る」の距離評価)が要る。
   ★`sorry` は使わず、**仮説として明示的に引数に取る**形にした
   (`Found/PGC/CohomologyColimit.lean` と同じ流儀)。
   ★★**偽ではない**: `PAdicLocalField p` は `Skeleton/PGC/Setup.lean:40` の定義により
   「`ℚ_[p]` の有限次拡大」であり、これは古典的な Ax–Sen–Tate の仮定
   (完備離散付値体・完全剰余体)を満たす。★つまり `AxSenTate K` は
   **真だが未証明**の命題であって、`hodgeTateDim_le_finrank` などが
   空虚に成り立っているのではない。
2. `Section3RealParameters.lean` の逸脱 3(Tate 捻りではなく固有空間で書いた)を
   そのまま引き継いでいる。
3. `hodgeTateDim_le_finrank` は `V` の**連続性を仮定していない**
   (`Representation` は連続性を持たない)。原典は連続表現に限る。
   ★不等号 `≤` の向きにはこの仮定は不要である。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped NormedField Valued Classical TensorProduct

/-! ## §1 抽象核

★★以下の 5 本には**分岐・付値・Galois・p 進の語彙が 1 つも出てこない**。 -/

section AbstractCore

variable {F A M ι G : Type*} [Field F] [Field A] [Algebra F A]
  [AddCommGroup M] [Module F M] [Module A M] [IsScalarTower F A M]

/-- ★★**抽象核** —— 半線型作用の固有ベクトルについての一次独立性の降下。

`A` に環準同型の族 `aut : G → A →+* A` が、`M` に加法準同型の族 `act : G → M →+ M` が
半線型(`hsem`)に作用し、

* `hfix` —— `aut` の全部で固定される `A` の元は `F` から来る、
* `hw` —— `w i` はどれも**同じ**固有値 `u g`(`≠ 0`)の固有ベクトル

とする。このとき `w` が `F`-一次独立なら `A`-一次独立である。

★`G` は**ただの型**でよい(群構造も `act` の乗法性も使わない)。
★`aut g` の単射性も使わない。 -/
theorem linearIndependent_of_fixedRange
    (act : G → M →+ M) (aut : G → A →+* A)
    (hsem : ∀ (g : G) (a : A) (m : M), act g (a • m) = aut g a • act g m)
    (hfix : ∀ a : A, (∀ g : G, aut g a = a) → a ∈ (algebraMap F A).range)
    {w : ι → M} {u : G → A} (hu : ∀ g : G, u g ≠ 0)
    (hw : ∀ (g : G) (i : ι), act g (w i) = u g • w i)
    (hFind : LinearIndependent F w) :
    LinearIndependent A w := by
  classical
  rw [linearIndependent_iff''] at hFind ⊢
  suffices H : ∀ n : ℕ, ∀ (s : Finset ι) (g : ι → A), s.card ≤ n →
      (∀ i ∉ s, g i = 0) → (∑ i ∈ s, g i • w i = 0) → ∀ i, g i = 0 by
    intro s g hs h0; exact H s.card s g le_rfl hs h0
  intro n
  induction n with
  | zero =>
      intro s g hcard hzero _ i
      have hs : s = ∅ := Finset.card_eq_zero.mp (Nat.le_zero.mp hcard)
      exact hzero i (by simp [hs])
  | succ n ih =>
      intro s g hcard hzero hsum i0
      by_contra hne
      have hi0s : i0 ∈ s := by
        by_contra hmem; exact hne (hzero i0 hmem)
      set c : ι → A := fun i => (g i0)⁻¹ * g i with hc
      have hci0 : c i0 = 1 := by simp [hc, inv_mul_cancel₀ hne]
      have hczero : ∀ i ∉ s, c i = 0 := fun i hi => by simp [hc, hzero i hi]
      have hcsum : ∑ i ∈ s, c i • w i = 0 := by
        have h2 : ∑ i ∈ s, c i • w i = (g i0)⁻¹ • ∑ i ∈ s, g i • w i := by
          rw [Finset.smul_sum]
          exact Finset.sum_congr rfl fun i _ => by rw [hc]; simp [mul_smul]
        rw [h2, hsum, smul_zero]
      -- 正規化した係数 `c` は `aut` の全部で固定される(台の最小性)
      have hfixed : ∀ (σ : G) (i : ι), aut σ (c i) = c i := by
        intro σ
        set d : ι → A := fun i => aut σ (c i) - c i with hd
        have hdi0 : d i0 = 0 := by simp [hd, hci0]
        have hdzero : ∀ i ∉ s.erase i0, d i = 0 := by
          intro i hi
          by_cases h : i = i0
          · subst h; exact hdi0
          · have hns : i ∉ s := fun hs => hi (Finset.mem_erase.mpr ⟨h, hs⟩)
            simp [hd, hczero i hns]
        have hautsum : ∑ i ∈ s, aut σ (c i) • w i = 0 := by
          have h1 : act σ (∑ i ∈ s, c i • w i) = u σ • ∑ i ∈ s, aut σ (c i) • w i := by
            rw [map_sum, Finset.smul_sum]
            exact Finset.sum_congr rfl fun i _ => by
              rw [hsem, hw, smul_smul, smul_smul, mul_comm]
          rw [hcsum, map_zero] at h1
          have h3 := congrArg (fun z => (u σ)⁻¹ • z) h1.symm
          simpa [smul_smul, inv_mul_cancel₀ (hu σ)] using h3
        have hdsum : ∑ i ∈ s.erase i0, d i • w i = 0 := by
          rw [Finset.sum_erase (f := fun i => d i • w i) s (by rw [hdi0, zero_smul])]
          have h4 : ∑ i ∈ s, d i • w i
              = (∑ i ∈ s, aut σ (c i) • w i) - ∑ i ∈ s, c i • w i := by
            rw [← Finset.sum_sub_distrib]
            exact Finset.sum_congr rfl fun i _ => by rw [hd]; exact sub_smul _ _ _
          rw [h4, hautsum, hcsum, sub_zero]
        have hcard' : (s.erase i0).card ≤ n := by
          have h5 := Finset.card_erase_of_mem hi0s
          omega
        intro i
        have h6 := ih (s.erase i0) d hcard' hdzero hdsum i
        simpa [hd, sub_eq_zero] using h6
      -- 固定される ⇒ `F` から来る ⇒ `F`-一次独立性に矛盾
      have hrange : ∀ i, ∃ b : F, algebraMap F A b = c i :=
        fun i => hfix (c i) (fun σ => hfixed σ i)
      choose b hb using hrange
      have hinj : Function.Injective (algebraMap F A) := (algebraMap F A).injective
      have hbzero : ∀ i ∉ s, b i = 0 := by
        intro i hi
        apply hinj
        rw [hb i, hczero i hi, map_zero]
      have hbsum : ∑ i ∈ s, b i • w i = 0 := by
        rw [← hcsum]
        exact Finset.sum_congr rfl fun i _ => by rw [← hb i, algebraMap_smul]
      have h7 := hFind s b hbzero hbsum i0
      rw [← hb i0, h7, map_zero] at hci0
      exact zero_ne_one hci0

/-- ★**抽象核の系** —— 固有ベクトルだけからなる `F`-部分空間の階数は
`M` の `A`-次元で押さえられる。 -/
theorem rank_le_finrank_of_fixedRange [Module.Finite A M]
    (act : G → M →+ M) (aut : G → A →+* A)
    (hsem : ∀ (g : G) (a : A) (m : M), act g (a • m) = aut g a • act g m)
    (hfix : ∀ a : A, (∀ g : G, aut g a = a) → a ∈ (algebraMap F A).range)
    {u : G → A} (hu : ∀ g : G, u g ≠ 0)
    (W : Submodule F M) (hW : ∀ x ∈ W, ∀ g : G, act g x = u g • x) :
    Module.rank F ↥W ≤ (Module.finrank A M : Cardinal) := by
  classical
  refine rank_le fun s hs => ?_
  have h1 : LinearIndependent F (fun i : ↥s => ((i : ↥W) : M)) :=
    hs.map' W.subtype (Submodule.ker_subtype W)
  have h2 : LinearIndependent A (fun i : ↥s => ((i : ↥W) : M)) :=
    linearIndependent_of_fixedRange act aut hsem hfix hu
      (fun g i => hW _ (i : ↥W).2 g) h1
  simpa using h2.fintype_card_le_finrank

/-- ★`finrank` 版。 -/
theorem finrank_le_of_fixedRange [Module.Finite A M]
    (act : G → M →+ M) (aut : G → A →+* A)
    (hsem : ∀ (g : G) (a : A) (m : M), act g (a • m) = aut g a • act g m)
    (hfix : ∀ a : A, (∀ g : G, aut g a = a) → a ∈ (algebraMap F A).range)
    {u : G → A} (hu : ∀ g : G, u g ≠ 0)
    (W : Submodule F M) (hW : ∀ x ∈ W, ∀ g : G, act g x = u g • x) :
    Module.finrank F ↥W ≤ Module.finrank A M :=
  Module.finrank_le_of_rank_le
    (rank_le_finrank_of_fixedRange act aut hsem hfix hu W hW)

/-- ★有限次元性。 -/
theorem finite_of_fixedRange [Module.Finite A M]
    (act : G → M →+ M) (aut : G → A →+* A)
    (hsem : ∀ (g : G) (a : A) (m : M), act g (a • m) = aut g a • act g m)
    (hfix : ∀ a : A, (∀ g : G, aut g a = a) → a ∈ (algebraMap F A).range)
    {u : G → A} (hu : ∀ g : G, u g ≠ 0)
    (W : Submodule F M) (hW : ∀ x ∈ W, ∀ g : G, act g x = u g • x) :
    Module.Finite F ↥W :=
  Module.rank_lt_aleph0_iff.mp
    (lt_of_le_of_lt (rank_le_finrank_of_fixedRange act aut hsem hfix hu W hW)
      Cardinal.natCast_lt_aleph0)

end AbstractCore

/-! ### 抽象核の非空虚性

★仮説が空虚でないことの証拠 —— `F = ℝ`, `A = ℂ`, `aut = 複素共役` で
上の仮説は**全部真になる**。 -/

/-- ★★**非空虚性** —— 実ベクトルの族が `ℝ` 上一次独立なら `ℂ` 上も一次独立。

抽象核 `linearIndependent_of_fixedRange` を `F = ℝ`, `A = ℂ`, `M = κ → ℂ`,
`G = Unit`(共役 1 個)、`u = 1` で使ったもの。
★これが通ること自体が「抽象核の仮説は空虚でない」ことの証拠である。 -/
theorem linearIndependent_complex_of_real {ι κ : Type*} {v : ι → κ → ℝ}
    (h : LinearIndependent ℝ fun i => (fun j => ((v i j : ℝ) : ℂ))) :
    LinearIndependent ℂ fun i => (fun j => ((v i j : ℝ) : ℂ)) := by
  refine linearIndependent_of_fixedRange (F := ℝ) (A := ℂ) (G := Unit)
    (fun _ => { toFun := fun m j => (starRingEnd ℂ) (m j)
                map_zero' := by ext j; simp
                map_add' := fun x y => by ext j; simp })
    (fun _ => starRingEnd ℂ) (fun _ a m => by ext j; simp) (fun a ha => ?_)
    (u := fun _ => (1 : ℂ)) (fun _ => one_ne_zero) (fun _ i => by ext j; simp) h
  exact ⟨a.re, by simpa using Complex.conj_eq_iff_re.mp (ha ())⟩

/-! ## §2 具体層 —— `ℂ_K^{Γ_K}` -/

variable {p : ℕ} [Fact p.Prime]

/-- `ℂ_K^{Γ_K}` —— `Γ_K` で固定される `ℂ_K` の元のなす `K`-部分空間。 -/
def fixedSubmodule (K : PAdicLocalField p) : Submodule K.carrier (CompKbar K) where
  carrier := {z | ∀ σ : K.absGal, σ • z = z}
  add_mem' {a b} ha hb σ := by rw [smul_add, ha σ, hb σ]
  zero_mem' σ := smul_zero σ
  smul_mem' c z hz σ := by rw [Algebra.smul_def, smul_mul', smul_algebraMap, hz σ]

@[simp] theorem mem_fixedSubmodule (K : PAdicLocalField p) (z : CompKbar K) :
    z ∈ fixedSubmodule K ↔ ∀ σ : K.absGal, σ • z = z := Iff.rfl

/-- ★**無条件に真** —— `K ⊆ ℂ_K^{Γ_K}`。Ax–Sen–Tate の易しい向き。 -/
theorem algebraMap_mem_fixedSubmodule (K : PAdicLocalField p) (c : K.carrier) :
    algebraMap K.carrier (CompKbar K) c ∈ fixedSubmodule K :=
  fun _ => smul_algebraMap _ _

/-- `1 ∈ ℂ_K^{Γ_K}`。 -/
theorem one_mem_fixedSubmodule (K : PAdicLocalField p) : (1 : CompKbar K) ∈ fixedSubmodule K :=
  fun σ => smul_one σ

/-- ★**非空虚性(具体層)** —— `ℂ_K^{Γ_K}` は `0` でない。

★`AxSenTate` は「空集合について何か言う」形の仮説ではない。 -/
theorem fixedSubmodule_ne_bot (K : PAdicLocalField p) : fixedSubmodule K ≠ ⊥ := by
  intro h
  exact one_ne_zero ((Submodule.mem_bot K.carrier).mp (h ▸ one_mem_fixedSubmodule K))

/-- ★★**Ax–Sen–Tate** —— `ℂ_K^{Γ_K} ⊆ K`。

★これは**仮説**であって、本ファイルは証明していない(逸脱の記録 1)。
逆向き `K ⊆ ℂ_K^{Γ_K}` は `algebraMap_mem_fixedSubmodule` で**無条件に真**なので、
この 1 本で `ℂ_K^{Γ_K} = K` と同値である。 -/
def AxSenTate (K : PAdicLocalField p) : Prop :=
  ∀ z : CompKbar K, (∀ σ : K.absGal, σ • z = z) →
    z ∈ (algebraMap K.carrier (CompKbar K)).range

def AxSenTate.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★`AxSenTate K` ⟺ `ℂ_K^{Γ_K}` の `K`-次元が `1`。 -/
theorem axSenTate_iff_finrank_fixedSubmodule (K : PAdicLocalField p) :
    AxSenTate K ↔ Module.finrank K.carrier ↥(fixedSubmodule K) = 1 := by
  have hone : (⟨1, one_mem_fixedSubmodule K⟩ : ↥(fixedSubmodule K)) ≠ 0 := by
    intro h
    exact one_ne_zero (congrArg Subtype.val h)
  rw [finrank_eq_one_iff_of_nonzero' _ hone]
  constructor
  · intro h z
    obtain ⟨c, hc⟩ := h (z : CompKbar K) z.2
    exact ⟨c, Subtype.ext (by simpa [Algebra.smul_def] using hc)⟩
  · intro h z hz
    obtain ⟨c, hc⟩ := h ⟨z, hz⟩
    exact ⟨c, by simpa [Algebra.smul_def] using congrArg Subtype.val hc⟩

/-! ## §3 `d_triv(0) = 1 ↔ Ax–Sen–Tate` -/

/-- `ℂ_K ⊗_{ℚ_p} ℚ_p ≃ ℂ_K` を `K`-線型同型として。 -/
noncomputable def tensorRid (K : PAdicLocalField p) :
    CompKbar K ⊗[ℚ_[p]] ℚ_[p] ≃ₗ[K.carrier] CompKbar K :=
  TensorProduct.AlgebraTensorModule.rid ℚ_[p] K.carrier (CompKbar K)

@[simp] theorem tensorRid_symm_apply (K : PAdicLocalField p) (z : CompKbar K) :
    (tensorRid K).symm z = z ⊗ₜ[ℚ_[p]] (1 : ℚ_[p]) :=
  TensorProduct.AlgebraTensorModule.rid_symm_apply _ _ _

/-- ★自明表現の重み `0` の空間は `tensorRid` で `ℂ_K^{Γ_K}` に移る。 -/
theorem map_weightSpace_zero (K : PAdicLocalField p) :
    Submodule.map (tensorRid K : _ →ₗ[K.carrier] CompKbar K)
        (hodgeTateWeightSpace K ℚ_[p] (Representation.trivial ℚ_[p] K.absGal ℚ_[p]) 0)
      = fixedSubmodule K := by
  ext z
  rw [Submodule.mem_map_equiv, mem_hodgeTateWeightSpace, mem_fixedSubmodule,
    tensorRid_symm_apply]
  constructor
  · intro h σ
    have := h σ
    simp only [Representation.tprod_apply, TensorProduct.map_tmul, compRep_apply,
      Representation.trivial_apply, cycloScalar_zero, one_smul] at this
    exact (tensorRid K).symm.injective (by
      rw [tensorRid_symm_apply, tensorRid_symm_apply]; exact this)
  · intro h σ
    simp only [Representation.tprod_apply, TensorProduct.map_tmul, compRep_apply,
      Representation.trivial_apply, cycloScalar_zero, one_smul, h σ]

/-- ★`d_triv(0)` は `ℂ_K^{Γ_K}` の `K`-次元そのもの(無条件)。 -/
theorem hodgeTateDim_trivial_zero (K : PAdicLocalField p) :
    hodgeTateDim K ℚ_[p] (Representation.trivial ℚ_[p] K.absGal ℚ_[p]) 0
      = Module.finrank K.carrier ↥(fixedSubmodule K) := by
  rw [hodgeTateDim,
    ((tensorRid K).submoduleMap
      (hodgeTateWeightSpace K ℚ_[p] (Representation.trivial ℚ_[p] K.absGal ℚ_[p]) 0)).finrank_eq]
  exact congrArg (fun S : Submodule K.carrier (CompKbar K) => Module.finrank K.carrier ↥S)
    (map_weightSpace_zero K)

/-- ★★★**本ファイルの主結果** —— `d_triv(0) = 1 ⟺ Ax–Sen–Tate`。

`Section3RealParameters.lean` が「`d_triv(0) = 1` すら言えない」と書いた原因は
**ちょうど Ax–Sen–Tate ただ 1 つ**であって、それ以外に足りないものは無い。
★`↔` なので、仮説が**弱すぎも強すぎもしない**ことまで言えている。 -/
theorem hodgeTateDim_trivial_zero_eq_one_iff (K : PAdicLocalField p) :
    hodgeTateDim K ℚ_[p] (Representation.trivial ℚ_[p] K.absGal ℚ_[p]) 0 = 1
      ↔ AxSenTate K := by
  rw [hodgeTateDim_trivial_zero, ← axSenTate_iff_finrank_fixedSubmodule]

def hodgeTateDim_trivial_zero_eq_one_iff.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-! ## §4 重み空間の有限次元性 -/

variable (K : PAdicLocalField p) (V : Type) [AddCommGroup V] [Module ℚ_[p] V]

/-- `Γ_K` の対角作用は `ℂ_K`-スカラーに対して**半線型**。 -/
theorem tprod_smul_semilinear (ρ : Representation ℚ_[p] K.absGal V) (σ : K.absGal)
    (a : CompKbar K) (z : CompKbar K ⊗[ℚ_[p]] V) :
    (Representation.tprod (compRep K) ρ) σ (a • z)
      = (σ • a) • (Representation.tprod (compRep K) ρ) σ z := by
  induction z using TensorProduct.induction_on with
  | zero => simp
  | tmul x v =>
      simp only [TensorProduct.smul_tmul', Representation.tprod_apply,
        TensorProduct.map_tmul, compRep_apply, smul_eq_mul, smul_mul']
  | add x y hx hy => rw [smul_add, map_add, map_add, hx, hy, smul_add]

/-- `χ_K(σ)^i` は `0` でない(`ℤ_[p]` の単数だから)。 -/
theorem cycloScalar_ne_zero (σ : K.absGal) (i : ℤ) : cycloScalar K σ i ≠ 0 := by
  rw [cycloScalar]
  simp

/-- ★★**`AxSenTate` を仮定すれば `d_V(i) ≤ dim_{ℚ_p} V`。**

★`Module.finrank` は無限次元で `0` を返すので、この不等式だけでは
「`d_V(i)` が意味を持つ」ことにはならない —— 一緒に
`finiteDimensional_hodgeTateWeightSpace` を見ること。 -/
theorem hodgeTateDim_le_finrank [FiniteDimensional ℚ_[p] V] (hK : AxSenTate K)
    (ρ : Representation ℚ_[p] K.absGal V) (i : ℤ) :
    hodgeTateDim K V ρ i ≤ Module.finrank ℚ_[p] V := by
  have h := finrank_le_of_fixedRange (F := K.carrier) (A := CompKbar K)
    (M := CompKbar K ⊗[ℚ_[p]] V) (G := K.absGal)
    (fun σ => ((Representation.tprod (compRep K) ρ) σ).toAddMonoidHom)
    (fun σ => MulSemiringAction.toRingHom K.absGal (CompKbar K) σ)
    (fun σ a z => tprod_smul_semilinear K V ρ σ a z) hK
    (u := fun σ => algebraMap ℚ_[p] (CompKbar K) (cycloScalar K σ i))
    (fun σ => by
      simpa using fun h => cycloScalar_ne_zero K σ i
        ((algebraMap ℚ_[p] (CompKbar K)).injective (by simpa using h)))
    (hodgeTateWeightSpace K V ρ i)
    (fun z hz σ => by
      simpa [algebraMap_smul] using hz σ)
  rwa [Module.finrank_baseChange (R := CompKbar K) (S := ℚ_[p]) (M' := V)] at h

def hodgeTateDim_le_finrank.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★**`AxSenTate` を仮定すれば重み空間は有限次元。** -/
theorem finiteDimensional_hodgeTateWeightSpace [FiniteDimensional ℚ_[p] V] (hK : AxSenTate K)
    (ρ : Representation ℚ_[p] K.absGal V) (i : ℤ) :
    FiniteDimensional K.carrier ↥(hodgeTateWeightSpace K V ρ i) :=
  finite_of_fixedRange (F := K.carrier) (A := CompKbar K) (M := CompKbar K ⊗[ℚ_[p]] V)
    (G := K.absGal)
    (fun σ => ((Representation.tprod (compRep K) ρ) σ).toAddMonoidHom)
    (fun σ => MulSemiringAction.toRingHom K.absGal (CompKbar K) σ)
    (fun σ a z => tprod_smul_semilinear K V ρ σ a z) hK
    (u := fun σ => algebraMap ℚ_[p] (CompKbar K) (cycloScalar K σ i))
    (fun σ => by
      simpa using fun h => cycloScalar_ne_zero K σ i
        ((algebraMap ℚ_[p] (CompKbar K)).injective (by simpa using h)))
    (hodgeTateWeightSpace K V ρ i)
    (fun z hz σ => by simpa [algebraMap_smul] using hz σ)

def finiteDimensional_hodgeTateWeightSpace.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★`AxSenTate` を仮定すると `d_triv(0) = 1`(主結果の片方向を明示したもの)。 -/
theorem hodgeTateDim_trivial_zero_eq_one (hK : AxSenTate K) :
    hodgeTateDim K ℚ_[p] (Representation.trivial ℚ_[p] K.absGal ℚ_[p]) 0 = 1 :=
  (hodgeTateDim_trivial_zero_eq_one_iff K).mpr hK

/-- ★★`IsHodgeTate` が**実際に満たされる例** —— `AxSenTate K` の下で
自明表現は Hodge-Tate(重みは `{0}` のみ)。

★`Section3RealParameters.lean` は `IsHodgeTate` を**定義しただけ**で、
それを満たす表現を 1 つも出せていなかった。★これが最初の 1 つである。 -/
theorem isHodgeTate_trivial (hK : AxSenTate K) :
    IsHodgeTate K ℚ_[p] (Representation.trivial ℚ_[p] K.absGal ℚ_[p]) :=
  ⟨{0}, by
    rw [Finset.sum_singleton, hodgeTateDim_trivial_zero_eq_one K hK, Module.finrank_self]⟩

def isHodgeTate_trivial.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-! ## §5 使っている公理の一覧

★どれも `propext` / `Classical.choice` / `Quot.sound` だけで、`sorryAx` は無い。
★`Classical.choice` が要るのは抽象核の `choose`(固定された係数から `F` の元を取る)と
`classical`(`DecidableEq ι`)である。★正規性・可判定性の追加公理は使っていない。 -/

#print axioms linearIndependent_of_fixedRange
#print axioms rank_le_finrank_of_fixedRange
#print axioms finrank_le_of_fixedRange
#print axioms finite_of_fixedRange
#print axioms linearIndependent_complex_of_real
#print axioms fixedSubmodule
#print axioms mem_fixedSubmodule
#print axioms algebraMap_mem_fixedSubmodule
#print axioms one_mem_fixedSubmodule
#print axioms AxSenTate
#print axioms axSenTate_iff_finrank_fixedSubmodule
#print axioms tensorRid
#print axioms tensorRid_symm_apply
#print axioms map_weightSpace_zero
#print axioms hodgeTateDim_trivial_zero
#print axioms hodgeTateDim_trivial_zero_eq_one_iff
#print axioms tprod_smul_semilinear
#print axioms cycloScalar_ne_zero
#print axioms hodgeTateDim_le_finrank
#print axioms finiteDimensional_hodgeTateWeightSpace
#print axioms hodgeTateDim_trivial_zero_eq_one
#print axioms fixedSubmodule_ne_bot
#print axioms isHodgeTate_trivial

end ABC3.Found.PGC
