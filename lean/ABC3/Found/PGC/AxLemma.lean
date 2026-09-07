import ABC3.Found.PGC.AxSenTate
import Mathlib.FieldTheory.Minpoly.IsConjRoot
import Mathlib.FieldTheory.Galois.Infinite

/-!
# [pGC] Ax の補題 —— `AxSenTate` を「代数元だけの距離評価」に落とす

`Found/PGC/AxSenTate.lean` は `d_triv(0) = 1 ↔ AxSenTate K` を証明し、
§3 の全部が **`AxSenTate K` ただ 1 つ**に依存することを `↔` で示した。
本ファイルはその `AxSenTate K` を **`ℂ_K`(完備化)の話から `K̄`(代数元)の話へ**
落とす。落とした先が古典的に **Ax の補題**と呼ばれるものである。

## 在庫調査(自分で測った。コマンドを残す)

```
grep -i -E "AxSen|Ax_Sen|SenTate|axSen" .cache/mathlib-index.txt      → 0 件
grep -i "hodgeTate" .cache/mathlib-index.txt                          → 0 件
grep -n "IsConjRoot" .cache/mathlib-index.txt                         → 40 件(★在る)
grep -n "fixedField" .cache/mathlib-index.txt                         → InfiniteGalois.fixedField_bot が在る
grep -n "isComplete_range" .cache/mathlib-index.txt                   → IsUniformInducing.isComplete_range が在る
grep -n "nextCoeff" .cache/mathlib-index.txt | grep -i "root"         → Splits.nextCoeff_eq_neg_sum_roots_of_monic
grep -n "	padicNorm\." .cache/mathlib-index.txt | grep -i one        → padicNorm.nat_eq_one_iff
```

★**Ax–Sen–Tate そのものは mathlib にも木にも無い**(直前の波の測定を追試して確認した)。
★一方、**その部品はほとんど mathlib に在った** —— 直前の波が「無い」と書いた語
(`AxSen`)ではなく、**型と名前空間**で引くと出る(`lean-idioms.md` #117(ii))。
とくに `IsConjRoot.exists_algEquiv`(`FieldTheory/Minpoly/IsConjRoot.lean:187`)は
「`minpoly` の根はちょうど `Γ_K` の軌道」を 1 行で与える。
★`import` を足さないと `Unknown identifier` になる(#68)——実測でそうなった。

## 何が言えたか

### §1 抽象核(分岐・付値・Galois・p 進の語彙が **1 つも出てこない**)

| 宣言 | 内容 |
|---|---|
| `norm_multisetSum_le_of_forall_le` | 超距離ノルム群の Multiset 和は各項の一様上界で押さえられる |
| `multisetSum_map_const_sub` | `Σ(x - r) = card • x - Σ r`(ただの加法群) |
| `norm_sub_multisetSum_div_le` | ★★**「重心」による近似**。誤差は `ε / ‖card‖` |
| `norm_sub_multisetSum_div_le_of_norm_card_eq_one` | ★ その系(`‖card‖ = 1` なら誤差は増えない) |
| `mem_of_forall_act_eq_of_approx` | ★★**「ほとんど不変ならほとんど `S`」⇒「不変なら `S`」** |

★`mem_of_forall_act_eq_of_approx` の `G` は**ただの型**でよい
(`AxSenTate.lean` の抽象核と同じく、群構造も作用の乗法性も使わない)。
★使うのは「`act g` が加法的」「`‖act g m‖ ≤ ‖m‖`」の 2 つだけである。
★`S` が部分群であることも要らない(**ただの閉集合**でよい)。

### §2 具体層

* `AxLemma K C` —— ★**Ax の補題**:
  `∀ σ ∈ Γ_K, ‖σx − x‖ ≤ ε` なる `x ∈ K̄` に対し `∃ y ∈ K, ‖x − y‖ ≤ C ε`。
* ★★`axSenTate_of_axLemma` —— **`AxLemma K C` (C ≥ 0) ⇒ `AxSenTate K`**。
  ★これが本ファイルの主結果。`AxSenTate` は `ℂ_K`(完備化)についての主張だが、
  それが**代数元だけの距離評価**に落ちる。
* ★`isClosed_range_algebraMap_compKbar` —— `K` は `ℂ_K` の中で閉
  (`K` が完備で `K → ℂ_K` が等長だから)。★これが「落とせる」ことの本体。

### §3 `AxLemma` のうち**証明できた部分**

★`AxLemma` 全体は証明していない(下記「残った穴」)。次は**証明した**:

0. ★★★`exists_norm_sub_algebraMap_le_div_norm_natDegree` ——
   **すべての `x` について、無条件に**

   `‖x − y‖ ≤ ε / ‖([K(x):K] : K)‖` なる `y ∈ K` が在る。

   ★これが本ファイルで到達した Ax の補題の一般形である。
   ★証明は「共役の重心 `(1/n) Σ_i x_i` を取る」だけで、
   **跡写像も Galois コホモロジーも分岐理論も 1 行も使わない**。
   ★★**足りないのは「定数が `x` に依らない」ことだけ**である。
1. ★★`exists_norm_sub_algebraMap_le_of_natDegree_tame` ——
   その系。**`[K(x):K]` が `p` と素なら `C = 1`**(tame な場合)。
   ★`exists_norm_sub_algebraMap_le_of_norm_natDegree_ge` はさらに一般で、
   `p^{v_p([K(x):K])} ≤ C` なる `x` について定数 `C` の Ax の補題を与える。
2. ★`axLemma_of_eps_zero` —— **`ε = 0` では無条件に成り立つ**(すべての `x` について)。
   実体は `mem_range_of_forall_smul_eq`(`K̄^{Γ_K} = K`)で、
   これは `InfiniteGalois.fixedField_bot` から出る。

★**非空虚性**: `natDegree_minpoly_algebraMap_tame` —— `x ∈ K` なら `[K(x):K] = 1` で
tame の仮定は満たされる。★つまり 1 の仮定は空集合について語っていない。

## ★残った穴(正直に書く)

★**`AxLemma K C` は埋まっていない。** 残っているのは
★★**「定数の一様性」ただ 1 点**である —— 上の 0 で得た定数 `‖([K(x):K] : K)‖⁻¹`
`= p^{v_p([K(x):K])}` は `x` について**非有界**なので、そのままでは
`axSenTate_of_axLemma` に入らない(近似列 `x_k` の次数が上がりうる)。

古典的には

* 次数 `p` の巡回拡大に落として、**共役差積(different)の評価**を使う(Sen の議論)、
* そこから `x` に依らない `C = |p|^{-p/(p-1)^2}` を得る(Ax)、

という段取りになる。★これは**数学が足りない**のであって配管ではない。
★必要な新ノード:「巡回 `p` 次拡大の different の評価(Sen の補題)」。

★★**`AxSenTate K` の項は作れていない。** 作れたのは
**`AxLemma K C → AxSenTate K`** という含意と、上の切片群である。

## 逸脱の記録(CLAUDE.md「逸脱」)

1. `AxLemma` の `ε` は **`0 ≤ ε`**(`0 < ε` ではない)で量化した。
   ★`ε = 0` を含めておくと `axLemma_of_eps_zero` がそのまま切片になる。
   原典(Ax)は `ε > 0` で書くが、`ε = 0` の場合は主張が強くなるだけで、
   後続(`axSenTate_of_axLemma`)は `ε > 0` しか使わないので影響しない。
2. tame の場合の仮定を `¬ p ∣ [K(x):K]` ではなく
   **`‖([K(x):K] : K)‖ = 1`** の形でも与えた(`..._of_norm_natDegree_eq_one`)。
   ★2 つは同値(`natDegree_tame_iff_norm_eq_one`)。証明で使うのは後者だけである。
3. `AxLemma` の `C` は「ある定数」として**外から与える**。原典は
   `C = |p|^{-p/(p-1)^2}` と具体値を出すが、`axSenTate_of_axLemma` は
   `0 ≤ C` しか使わない。★値を決めないほうが仮説として弱い。
4. 本ファイルの `.src` は `AxSenTate.lean` と**同じ項目**(pGC 物理 p.6 Corollary 3.1)を
   指している。本ファイルの中身は `AxSenTate` の**入力**であって、原典が独立の項目として
   立てているものではないため。★原典に無い分解であることを明示しておく。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped NormedField Valued

/-! ## §1 抽象核

★★以下の 5 本には**分岐・付値・Galois・p 進の語彙が 1 つも出てこない**。 -/

section AbstractCore

/-- **抽象核** —— 超距離ノルム群では、Multiset 和のノルムは各項の一様上界で押さえられる。

★項数に依存しない(これがアルキメデス的な三角不等式との違い)。 -/
theorem norm_multisetSum_le_of_forall_le {M : Type*} [SeminormedAddCommGroup M]
    [IsUltrametricDist M] {s : Multiset M} {C : ℝ}
    (hC : 0 ≤ C) (h : ∀ a ∈ s, ‖a‖ ≤ C) : ‖s.sum‖ ≤ C := by
  induction s using Multiset.induction with
  | empty => simpa using hC
  | cons a t ih =>
      rw [Multiset.sum_cons]
      refine (IsUltrametricDist.norm_add_le_max a t.sum).trans (max_le ?_ ?_)
      · exact h a (Multiset.mem_cons_self a t)
      · exact ih fun b hb => h b (Multiset.mem_cons_of_mem hb)

/-- **抽象核** —— `Σ_{r ∈ s} (x - r) = card s • x - Σ_{r ∈ s} r`。ただの加法群でよい。 -/
theorem multisetSum_map_const_sub {A : Type*} [AddCommGroup A] (s : Multiset A) (x : A) :
    (s.map (fun r => x - r)).sum = Multiset.card s • x - s.sum := by
  induction s using Multiset.induction with
  | empty => simp
  | cons a t ih =>
      rw [Multiset.map_cons, Multiset.sum_cons, ih, Multiset.sum_cons, Multiset.card_cons,
        succ_nsmul]
      abel

/-- ★★★**抽象核** —— 「重心による近似」。

`s` の各元が `x` から `ε` 以内にあるとき、重心 `(Σ s)/card s` は
`x` から **`ε / ‖card s‖` 以内**にある。

★超距離性が効くのはここ 1 箇所である(**項数分の増幅が起きない**)。
★`card s` で割る分だけ誤差が `‖card s‖⁻¹` 倍に増える —— これが
Ax の補題の定数の出所であり、★**そこだけが残った穴**である
(§3 の docstring を見ること)。 -/
theorem norm_sub_multisetSum_div_le {A : Type*} [NormedField A] [IsUltrametricDist A]
    {s : Multiset A} {x : A} {ε : ℝ} (hε : 0 ≤ ε)
    (hne : (Multiset.card s : A) ≠ 0)
    (h : ∀ r ∈ s, ‖x - r‖ ≤ ε) :
    ‖x - s.sum / (Multiset.card s : A)‖ ≤ ε / ‖(Multiset.card s : A)‖ := by
  have key : x - s.sum / (Multiset.card s : A)
      = (s.map (fun r => x - r)).sum / (Multiset.card s : A) := by
    rw [multisetSum_map_const_sub, nsmul_eq_mul]
    field_simp
  have h2 : ‖(s.map (fun r => x - r)).sum‖ ≤ ε := by
    refine norm_multisetSum_le_of_forall_le hε ?_
    intro a ha
    obtain ⟨r, hr, rfl⟩ := Multiset.mem_map.mp ha
    exact h r hr
  rw [key, norm_div, div_eq_mul_inv, div_eq_mul_inv]
  exact mul_le_mul_of_nonneg_right h2 (by positivity)

/-- ★**抽象核の系** —— `‖card s‖ = 1`(「tame」)なら誤差は増えない。 -/
theorem norm_sub_multisetSum_div_le_of_norm_card_eq_one {A : Type*} [NormedField A]
    [IsUltrametricDist A] {s : Multiset A} {x : A} {ε : ℝ} (hε : 0 ≤ ε)
    (hn : ‖(Multiset.card s : A)‖ = 1)
    (h : ∀ r ∈ s, ‖x - r‖ ≤ ε) :
    ‖x - s.sum / (Multiset.card s : A)‖ ≤ ε := by
  have hne : (Multiset.card s : A) ≠ 0 := by
    intro h0; rw [h0, norm_zero] at hn; exact zero_ne_one hn
  have := norm_sub_multisetSum_div_le hε hne h
  rwa [hn, div_one] at this

/-- ★★★**抽象核** —— 「ほとんど不変な元がほとんど `S` に居る」なら
「不変な元は `S` に居る」。

* `act` —— 型 `G` で添字された**加法準同型**の族。ノルムを増やさない(`hact`)。
* `S` —— **閉集合**(部分群でなくてよい)。
* `D` —— **稠密集合**。
* `happrox` —— `D` の元については「ほとんど不変 ⇒ ほとんど `S`」が定数 `C` で成り立つ。

★`G` は**ただの型**でよい(群構造も `act` の乗法性も使わない)。
★超距離性も要らない(ふつうの三角不等式だけで通る)。 -/
theorem mem_of_forall_act_eq_of_approx {M : Type*} [SeminormedAddCommGroup M] {G : Type*}
    (act : G → M →+ M) (hact : ∀ (g : G) (m : M), ‖act g m‖ ≤ ‖m‖)
    {S : Set M} (hS : IsClosed S) {D : Set M} (hD : Dense D)
    {C : ℝ} (hC : 0 ≤ C)
    (happrox : ∀ x ∈ D, ∀ ε : ℝ, 0 ≤ ε → (∀ g : G, ‖act g x - x‖ ≤ ε) →
      ∃ y ∈ S, ‖x - y‖ ≤ C * ε)
    {z : M} (hz : ∀ g : G, act g z = z) : z ∈ S := by
  rw [← hS.closure_eq, Metric.mem_closure_iff]
  intro ε hε
  have hden : (0 : ℝ) < 2 * C + 2 := by linarith
  set δ : ℝ := ε / (2 * C + 2) with hδdef
  have hδ : 0 < δ := div_pos hε hden
  obtain ⟨x, hxD, hx⟩ := Metric.mem_closure_iff.mp (hD z) δ hδ
  have hzx : ‖z - x‖ < δ := by rwa [← dist_eq_norm]
  have hstep : ∀ g : G, ‖act g x - x‖ ≤ 2 * δ := by
    intro g
    have h1 : act g x - x = act g (x - z) + (z - x) := by
      rw [map_sub, hz g]; abel
    have h2 : ‖act g (x - z)‖ ≤ ‖x - z‖ := hact g _
    have h3 : ‖x - z‖ = ‖z - x‖ := norm_sub_rev _ _
    have h4 := norm_add_le (act g (x - z)) (z - x)
    rw [← h1] at h4
    linarith
  obtain ⟨y, hyS, hy⟩ := happrox x hxD (2 * δ) (by positivity) hstep
  refine ⟨y, hyS, ?_⟩
  have htri : ‖z - y‖ ≤ ‖z - x‖ + ‖x - y‖ := by
    simpa using norm_add_le (z - x) (x - y)
  have h5 : C * (2 * δ) = 2 * (C * δ) := by ring
  rw [h5] at hy
  have hδeq : ε = 2 * (C * δ) + 2 * δ := by
    rw [hδdef]; field_simp
  rw [dist_eq_norm]
  linarith

end AbstractCore

/-! ## §2 具体層 —— `AxLemma K C ⇒ AxSenTate K` -/

variable {p : ℕ} [Fact p.Prime]

/-- ★★**Ax の補題**(距離評価)。

「すべての `σ ∈ Γ_K` について `‖σx − x‖ ≤ ε`」ならば
「`‖x − y‖ ≤ C ε` なる `y ∈ K` が在る」。

★これは**仮説**であって、本ファイルは全体としては証明していない。
証明できた切片は `exists_norm_sub_algebraMap_le_of_natDegree_tame`(tame な場合)と
`axLemma_of_eps_zero`(`ε = 0`)の 2 枚である。 -/
def AxLemma (K : PAdicLocalField p) (C : ℝ) : Prop :=
  ∀ ε : ℝ, 0 ≤ ε → ∀ x : K.closure, (∀ σ : K.absGal, ‖σ • x - x‖ ≤ ε) →
    ∃ y : K.carrier, ‖x - algebraMap K.carrier K.closure y‖ ≤ C * ε

def AxLemma.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- `Γ_K` の `ℂ_K` への作用を加法準同型として見たもの。 -/
noncomputable def smulAddHom (K : PAdicLocalField p) (σ : K.absGal) :
    CompKbar K →+ CompKbar K :=
  DistribSMul.toAddMonoidHom (CompKbar K) σ

@[simp] theorem smulAddHom_apply (K : PAdicLocalField p) (σ : K.absGal) (w : CompKbar K) :
    smulAddHom K σ w = σ • w := rfl

/-- `K ↪ ℂ_K` は等長。 -/
theorem isometry_algebraMap_compKbar (K : PAdicLocalField p) :
    Isometry (algebraMap K.carrier (CompKbar K)) :=
  Isometry.of_dist_eq fun a b => by
    rw [dist_eq_norm, dist_eq_norm, ← map_sub, norm_algebraMap_closureCompletion]

/-- ★★`K` の像は `ℂ_K` の中で**閉**。

`K` は `ℚ_[p]` 上有限次だから完備で、`K → ℂ_K` は等長なので像は完備、したがって閉。
★これが「`ℂ_K` の話を `K̄` の話に落とせる」ことの本体である。 -/
theorem isClosed_range_algebraMap_compKbar (K : PAdicLocalField p) :
    IsClosed (Set.range (algebraMap K.carrier (CompKbar K))) :=
  ((isometry_algebraMap_compKbar K).isUniformInducing.isComplete_range).isClosed

/-- ★★★**本ファイルの主結果** —— **Ax の補題から Ax–Sen–Tate が出る**。

`AxSenTate K` は完備化 `ℂ_K` についての主張だが、`K̄` が `ℂ_K` で稠密で
`K` が `ℂ_K` で閉であることから、**代数元だけの距離評価**(`AxLemma`)に落ちる。 -/
theorem axSenTate_of_axLemma {K : PAdicLocalField p} {C : ℝ} (hC : 0 ≤ C)
    (h : AxLemma K C) : AxSenTate K := by
  intro z hz
  have hmem : z ∈ Set.range (algebraMap K.carrier (CompKbar K)) := by
    refine mem_of_forall_act_eq_of_approx (smulAddHom K)
      (fun σ m => le_of_eq (norm_smul_closureCompletion K σ m))
      (isClosed_range_algebraMap_compKbar K)
      (denseRange_coe_closureCompletion K) hC ?_ (fun σ => hz σ)
    rintro _ ⟨x0, rfl⟩ ε hε hact
    have hx : ∀ σ : K.absGal, ‖σ • x0 - x0‖ ≤ ε := by
      intro σ
      have h1 := hact σ
      rw [smulAddHom_apply] at h1
      have h2 : σ • (x0 : CompKbar K) - (x0 : CompKbar K)
          = ((σ • x0 - x0 : K.closure) : CompKbar K) := by
        push_cast
        rw [smul_coe_closureCompletion]
      rw [h2, norm_coe_closureCompletion] at h1
      exact h1
    obtain ⟨y, hy⟩ := h ε hε x0 hx
    refine ⟨algebraMap K.carrier (CompKbar K) y, ⟨y, rfl⟩, ?_⟩
    have h3 : (x0 : CompKbar K) - algebraMap K.carrier (CompKbar K) y
        = ((x0 - algebraMap K.carrier K.closure y : K.closure) : CompKbar K) := by
      rw [UniformSpace.Completion.algebraMap_def]
      push_cast
      rfl
    rw [h3, norm_coe_closureCompletion]
    exact hy
  obtain ⟨c, hc⟩ := hmem
  exact ⟨c, hc⟩

def axSenTate_of_axLemma.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-! ## §3 `AxLemma` の切片 1 —— tame(`p ∤ [K(x):K]`)な場合 -/

/-- `‖(n : K̄)‖ = ‖(n : K)‖`。 -/
theorem norm_natCast_closure (K : PAdicLocalField p) (n : ℕ) :
    ‖((n : ℕ) : K.closure)‖ = ‖((n : ℕ) : K.carrier)‖ := by
  rw [show ((n : ℕ) : K.closure) = algebraMap K.carrier K.closure ((n : ℕ) : K.carrier) from
    (map_natCast _ n).symm, norm_algebraMap_closure]

/-- `‖(n : K)‖ = ‖(n : ℚ_p)‖`(スペクトルノルムは `ℚ_p` のノルムを延長する)。 -/
theorem norm_natCast_carrier (K : PAdicLocalField p) (n : ℕ) :
    ‖((n : ℕ) : K.carrier)‖ = ‖((n : ℕ) : ℚ_[p])‖ := by
  rw [show ((n : ℕ) : K.carrier) = algebraMap ℚ_[p] K.carrier ((n : ℕ) : ℚ_[p]) from
    (map_natCast _ n).symm, norm_algebraMap]

/-- `n ≠ 0` なら `‖(n : K̄)‖ > 0`(`K̄` は標数 0)。 -/
theorem norm_natCast_closure_pos (K : PAdicLocalField p) {n : ℕ} (hn : n ≠ 0) :
    0 < ‖((n : ℕ) : K.closure)‖ := by
  rw [norm_natCast_closure, norm_natCast_carrier]
  exact norm_pos_iff.mpr (Nat.cast_ne_zero.mpr hn)

/-- ★`‖(n : K)‖ = 1 ⟺ `p ∤ n`。「tame」の言い換え。 -/
theorem natDegree_tame_iff_norm_eq_one (K : PAdicLocalField p) (n : ℕ) :
    ‖((n : ℕ) : K.carrier)‖ = 1 ↔ ¬ (p ∣ n) := by
  rw [norm_natCast_carrier, ← padicNorm.nat_eq_one_iff (p := p) n,
    show ((n : ℕ) : ℚ_[p]) = (((n : ℕ) : ℚ) : ℚ_[p]) by push_cast; ring,
    ← Padic.padicNormE.is_norm, padicNormE.eq_padic_norm']
  exact_mod_cast Iff.rfl

/-- ★`x` の共役(`minpoly` の根)はすべて `x` から `ε` 以内。

`minpoly` の根はちょうど `Γ_K` の軌道である(`IsConjRoot.exists_algEquiv`)。 -/
theorem norm_sub_le_of_mem_aroots (K : PAdicLocalField p) {x r : K.closure} {ε : ℝ}
    (hx : ∀ σ : K.absGal, ‖σ • x - x‖ ≤ ε)
    (hr : r ∈ (minpoly K.carrier x).aroots K.closure) : ‖x - r‖ ≤ ε := by
  have hint : IsIntegral K.carrier x :=
    (Algebra.IsAlgebraic.isAlgebraic (R := K.carrier) x).isIntegral
  have hconj : IsConjRoot K.carrier x r := (isConjRoot_iff_mem_minpoly_aroots hint).mpr hr
  obtain ⟨σ, hσ⟩ := hconj.exists_algEquiv
  have hr' : σ⁻¹ • x = r := by
    rw [← hσ]
    show σ⁻¹ • (σ • r) = r
    rw [inv_smul_smul]
  rw [← hr', norm_sub_rev]
  exact hx σ⁻¹

/-- ★★★**Ax の補題の「証明できた一般形」** ——
`∀ σ ∈ Γ_K, ‖σx − x‖ ≤ ε` なら、`y ∈ K` を取って

  `‖x − y‖ ≤ ε / ‖([K(x):K] : K)‖`

にできる。★**定数が `x` の次数に依存している**のがこの版の限界である
(古典的な Ax の定理は、これを **`x` に依らない定数 `C`** に置き換える)。

証明は「共役の重心 `(1/n) Σ_i x_i` を取る」だけ:
★`Σ_i x_i` は `minpoly` の `nextCoeff` なので **`K` の元**であり、
共役はすべて `x` から `ε` 以内(`norm_sub_le_of_mem_aroots`)、
超距離性から和も `ε` 以内、`n` で割って `ε/‖n‖`。
★**跡写像も Galois コホモロジーも分岐理論も 1 行も使っていない。** -/
theorem exists_norm_sub_algebraMap_le_div_norm_natDegree
    (K : PAdicLocalField p) {x : K.closure} {ε : ℝ} (hε : 0 ≤ ε)
    (hx : ∀ σ : K.absGal, ‖σ • x - x‖ ≤ ε) :
    ∃ y : K.carrier, ‖x - algebraMap K.carrier K.closure y‖
      ≤ ε / ‖(((minpoly K.carrier x).natDegree : ℕ) : K.carrier)‖ := by
  classical
  have hint : IsIntegral K.carrier x :=
    (Algebra.IsAlgebraic.isAlgebraic (R := K.carrier) x).isIntegral
  have hmonic : (minpoly K.carrier x).Monic := minpoly.monic hint
  have hdeg : (minpoly K.carrier x).natDegree ≠ 0 := (minpoly.natDegree_pos hint).ne'
  set g : Polynomial K.closure :=
    (minpoly K.carrier x).map (algebraMap K.carrier K.closure) with hg
  have hgmonic : g.Monic := hmonic.map _
  have hsplits : g.Splits := IsAlgClosed.splits g
  have hcard : Multiset.card g.roots = (minpoly K.carrier x).natDegree := by
    rw [Polynomial.splits_iff_card_roots.mp hsplits, hg, hmonic.natDegree_map]
  have hnext : g.roots.sum
      = algebraMap K.carrier K.closure (-(minpoly K.carrier x).nextCoeff) := by
    have h1 := hsplits.nextCoeff_eq_neg_sum_roots_of_monic hgmonic
    rw [hg, Polynomial.nextCoeff_map_eq] at h1
    rw [map_neg, h1, neg_neg]
  have hne : ((Multiset.card g.roots : ℕ) : K.closure) ≠ 0 := by
    rw [hcard]
    exact norm_pos_iff.mp (norm_natCast_closure_pos K hdeg)
  have hmain := norm_sub_multisetSum_div_le (x := x) hε hne
    (fun r hr => norm_sub_le_of_mem_aroots K hx hr)
  refine ⟨(-(minpoly K.carrier x).nextCoeff) / (((minpoly K.carrier x).natDegree : ℕ) :
    K.carrier), ?_⟩
  have hmap : algebraMap K.carrier K.closure
        ((-(minpoly K.carrier x).nextCoeff) / (((minpoly K.carrier x).natDegree : ℕ) : K.carrier))
      = g.roots.sum / ((Multiset.card g.roots : ℕ) : K.closure) := by
    rw [map_div₀, hnext, hcard, map_natCast]
  rw [hmap, hcard]
  rwa [hcard, norm_natCast_closure] at hmain

def exists_norm_sub_algebraMap_le_div_norm_natDegree.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★**Ax の補題の tame な切片(ノルムの形)** —— `‖([K(x):K] : K)‖ = 1` なら `C = 1`。 -/
theorem exists_norm_sub_algebraMap_le_of_norm_natDegree_eq_one
    (K : PAdicLocalField p) {x : K.closure} {ε : ℝ} (hε : 0 ≤ ε)
    (hn : ‖(((minpoly K.carrier x).natDegree : ℕ) : K.carrier)‖ = 1)
    (hx : ∀ σ : K.absGal, ‖σ • x - x‖ ≤ ε) :
    ∃ y : K.carrier, ‖x - algebraMap K.carrier K.closure y‖ ≤ ε := by
  obtain ⟨y, hy⟩ := exists_norm_sub_algebraMap_le_div_norm_natDegree K hε hx
  exact ⟨y, by rwa [hn, div_one] at hy⟩

def exists_norm_sub_algebraMap_le_of_norm_natDegree_eq_one.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★**`AxLemma K C` は「次数の `p` 冪部分が `C` で抑えられる `x`」については成り立つ**。

`1 ≤ C * ‖([K(x):K] : K)‖`、すなわち `p^{v_p([K(x):K])} ≤ C` のとき、
その `x` については定数 `C` の Ax の補題が成り立つ。
★★**したがって残っている穴は「`v_p([K(x):K])` が非有界であること」ただ 1 点である。** -/
theorem exists_norm_sub_algebraMap_le_of_norm_natDegree_ge
    (K : PAdicLocalField p) {x : K.closure} {ε C : ℝ} (hε : 0 ≤ ε)
    (hC : 1 ≤ C * ‖(((minpoly K.carrier x).natDegree : ℕ) : K.carrier)‖)
    (hx : ∀ σ : K.absGal, ‖σ • x - x‖ ≤ ε) :
    ∃ y : K.carrier, ‖x - algebraMap K.carrier K.closure y‖ ≤ C * ε := by
  have hint : IsIntegral K.carrier x :=
    (Algebra.IsAlgebraic.isAlgebraic (R := K.carrier) x).isIntegral
  have hdeg : (minpoly K.carrier x).natDegree ≠ 0 := (minpoly.natDegree_pos hint).ne'
  have hnpos : 0 < ‖(((minpoly K.carrier x).natDegree : ℕ) : K.carrier)‖ := by
    rw [norm_natCast_carrier]
    exact norm_pos_iff.mpr (Nat.cast_ne_zero.mpr hdeg)
  obtain ⟨y, hy⟩ := exists_norm_sub_algebraMap_le_div_norm_natDegree K hε hx
  refine ⟨y, hy.trans ?_⟩
  have h1 : 1 / ‖(((minpoly K.carrier x).natDegree : ℕ) : K.carrier)‖ ≤ C := by
    rw [div_le_iff₀ hnpos]; exact hC
  calc ε / ‖(((minpoly K.carrier x).natDegree : ℕ) : K.carrier)‖
      = ε * (1 / ‖(((minpoly K.carrier x).natDegree : ℕ) : K.carrier)‖) := by ring
    _ ≤ ε * C := mul_le_mul_of_nonneg_left h1 hε
    _ = C * ε := mul_comm _ _

def exists_norm_sub_algebraMap_le_of_norm_natDegree_ge.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★**Ax の補題の tame な切片(可読な形)** —— `p ∤ [K(x):K]` なら `C = 1`。 -/
theorem exists_norm_sub_algebraMap_le_of_natDegree_tame
    (K : PAdicLocalField p) {x : K.closure} {ε : ℝ} (hε : 0 ≤ ε)
    (hn : ¬ (p ∣ (minpoly K.carrier x).natDegree))
    (hx : ∀ σ : K.absGal, ‖σ • x - x‖ ≤ ε) :
    ∃ y : K.carrier, ‖x - algebraMap K.carrier K.closure y‖ ≤ ε :=
  exists_norm_sub_algebraMap_le_of_norm_natDegree_eq_one K hε
    ((natDegree_tame_iff_norm_eq_one K _).mpr hn) hx

def exists_norm_sub_algebraMap_le_of_natDegree_tame.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★**非空虚性** —— `x ∈ K` なら `[K(x):K] = 1` なので tame の仮定は満たされる。

★`exists_norm_sub_algebraMap_le_of_natDegree_tame` は空集合について語っていない。 -/
theorem natDegree_minpoly_algebraMap_tame (K : PAdicLocalField p) (c : K.carrier) :
    ¬ (p ∣ (minpoly K.carrier (algebraMap K.carrier K.closure c)).natDegree) := by
  rw [minpoly.eq_X_sub_C K.closure c, Polynomial.natDegree_X_sub_C, Nat.dvd_one]
  exact (Fact.out : p.Prime).one_lt.ne'

/-! ## §4 `AxLemma` の切片 2 —— `ε = 0`(`K̄^{Γ_K} = K`) -/

/-- ★★**`K̄^{Γ_K} = K`** —— 完備化を取る前の(易しい方の)Ax–Sen–Tate。

`K̄/K` は Galois なので `InfiniteGalois.fixedField_bot` がそのまま効く。
★**完備化を取った `ℂ_K` では同じ議論が効かない** —— そこが Ax の補題の要る所である。 -/
theorem mem_range_of_forall_smul_eq (K : PAdicLocalField p) {x : K.closure}
    (hx : ∀ σ : K.absGal, σ • x = x) :
    x ∈ (algebraMap K.carrier K.closure).range := by
  have hmem : x ∈ IntermediateField.fixedField
      (⊤ : Subgroup (K.closure ≃ₐ[K.carrier] K.closure)) := by
    rw [IntermediateField.mem_fixedField_iff]
    intro σ _
    exact hx σ
  rw [InfiniteGalois.fixedField_bot] at hmem
  obtain ⟨c, hc⟩ := IntermediateField.mem_bot.mp hmem
  exact ⟨c, hc⟩

def mem_range_of_forall_smul_eq.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★**Ax の補題は `ε = 0` では無条件に成り立つ**(定数 `C` は何でもよい)。 -/
theorem axLemma_of_eps_zero (K : PAdicLocalField p) (C : ℝ) (x : K.closure)
    (hx : ∀ σ : K.absGal, ‖σ • x - x‖ ≤ 0) :
    ∃ y : K.carrier, ‖x - algebraMap K.carrier K.closure y‖ ≤ C * 0 := by
  have hfix : ∀ σ : K.absGal, σ • x = x := by
    intro σ
    have := hx σ
    have h0 : ‖σ • x - x‖ = 0 := le_antisymm this (norm_nonneg _)
    exact sub_eq_zero.mp (norm_eq_zero.mp h0)
  obtain ⟨y, hy⟩ := mem_range_of_forall_smul_eq K hfix
  exact ⟨y, by rw [hy, sub_self, norm_zero, mul_zero]⟩

def axLemma_of_eps_zero.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-! ## §5 使っている公理の一覧

★どれも `propext` / `Classical.choice` / `Quot.sound` だけで、`sorryAx` は**無い**。 -/

#print axioms norm_multisetSum_le_of_forall_le
#print axioms multisetSum_map_const_sub
#print axioms norm_sub_multisetSum_div_le
#print axioms mem_of_forall_act_eq_of_approx
#print axioms AxLemma
#print axioms smulAddHom
#print axioms smulAddHom_apply
#print axioms isometry_algebraMap_compKbar
#print axioms isClosed_range_algebraMap_compKbar
#print axioms axSenTate_of_axLemma
#print axioms norm_natCast_closure
#print axioms norm_natCast_carrier
#print axioms natDegree_tame_iff_norm_eq_one
#print axioms norm_sub_le_of_mem_aroots
#print axioms norm_sub_multisetSum_div_le_of_norm_card_eq_one
#print axioms norm_natCast_closure_pos
#print axioms exists_norm_sub_algebraMap_le_div_norm_natDegree
#print axioms exists_norm_sub_algebraMap_le_of_norm_natDegree_eq_one
#print axioms exists_norm_sub_algebraMap_le_of_norm_natDegree_ge
#print axioms exists_norm_sub_algebraMap_le_of_natDegree_tame
#print axioms natDegree_minpoly_algebraMap_tame
#print axioms mem_range_of_forall_smul_eq
#print axioms axLemma_of_eps_zero

end ABC3.Found.PGC
