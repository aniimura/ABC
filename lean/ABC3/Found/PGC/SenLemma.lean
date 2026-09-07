import ABC3.Found.PGC.AxLemma
import Mathlib.LinearAlgebra.Lagrange

/-!
# [pGC] Sen の補題に向けて —— 「重み付き重心」と「降下」

`Found/PGC/AxLemma.lean` は `AxLemma K C → AxSenTate K` を証明し、
`AxLemma K C` のうち **tame な切片**(`p ∤ [K(x):K]` なら `C = 1`)と
**無条件だが定数が `x` に依存する形**

```
exists_norm_sub_algebraMap_le_div_norm_natDegree :
  ∃ y ∈ K, ‖x − y‖ ≤ ε / ‖([K(x):K] : K)‖
```

まで到達した。残っていたのは **「定数の一様性」ただ 1 点**である。
本ファイルはそこを詰めるための道具を 2 種類そろえる。

## 在庫調査(自分で測った。コマンドを残す。★MCP は 0 回)

```
grep -n "IsUltrametricDist" .cache/mathlib-index.txt | grep -i "sum"      → 0 件
grep -n "norm_sum_le_of_forall_le\|nnorm_sum_le" .cache/mathlib-index.txt → 0 件
   ⇒ 超距離での Multiset 和の一様上界は mathlib に無い(AxLemma.lean が自前で持つ)
grep -n "Polynomial.roots_map\|roots_map_of\|aroots_map" .cache/mathlib-index.txt
   → ★Polynomial.Splits.roots_map_of_injective (Algebra/Polynomial/Splits.lean:392) が在る
grep -n "aeval_algHom_apply" .cache/mathlib-index.txt
   → ★Polynomial.aeval_algHom_apply (Algebra/Polynomial/AlgebraMap.lean:413) が在る
grep -n "minpoly.*natDegree_eq_one" .cache/mathlib-index.txt
   → ★minpoly.natDegree_eq_one_iff (FieldTheory/Minpoly/Basic.lean:220) が在る
grep -n "differentIdeal" .cache/mathlib-index.txt                        → 20 件(Dedekind 版)
grep -n "traceDual\|FractionalIdeal.dual" .cache/mathlib-index.txt       → ★在る
   ★Module.Basis.traceDual_powerBasis_eq (RingTheory/Trace/Basic.lean:610):
     pb.basis.traceDual i = (minpolyDiv K pb.gen).coeff i / aeval pb.gen (derivative (minpoly K pb.gen))
   ——Euler の等式 `Σ_r r^{n-1}/f'(r) = 1` の mathlib 側の入口だが `PowerBasis` を要求する。
grep -n "Lagrange\." .cache/mathlib-index.txt
   → ★★Lagrange.coeff_eq_sum (LinearAlgebra/Lagrange.lean:495):
       P.coeff (#s − 1) = Σ_{i ∈ s} P.eval (v i) / Π_{j ≠ i} (v i − v j)
     ★★**こちらの方が軽い** —— 中間体も PowerBasis も要らず、根の集合の上だけで済む。
     ★本ファイルは Euler の等式をこちらで証明した(`aroots_sum_div_derivative_eq_one`)。
grep -n "eval_multiset_prod_X_sub_C_derivative" .cache/mathlib-index.txt
   → ★在る(Algebra/Polynomial/Derivative.lean:680)。`Π_{j≠r}(r−j) = f′(r)` はこれ 1 行。
awk -F'\t' '$2 ~ /ABC3.*smul.*closure/' .cache/decl-index.txt
   → ★smul_closure_def / norm_smul_closure (Found/PGC/AbsClosureModules.lean:260,270)
```

★**「Ax–Sen–Tate が無い」は真だが「部品が無い」は偽**、という直前の波の教訓は
本ファイルでも再現した ——「根の重複集合が Galois 安定」は
**`Polynomial.Splits.roots_map_of_injective` の 1 行**で出る。

## 何が言えたか

### §1 抽象核(分岐・付値・Galois・p 進の語彙が **1 つも出てこない**)

| 宣言 | 内容 |
|---|---|
| `multisetSum_weighted_const_sub` | `Σ w(r)(x−r) = (Σ w(r))·x − Σ w(r)r`(ただの可換環) |
| `norm_sub_multisetSum_weighted_le` | ★★**重み付き重心**。`Σ w = 1`・`‖w‖ ≤ B` なら誤差は `B ε` |
| `exists_mem_of_descent` | ★★**超距離での降下**。`deg` が下がる各段の誤差が `ε` なら合計も `ε` |

★`exists_mem_of_descent` は `‖·‖` と `ℕ` しか使わない。強帰納法 1 本。
★超距離性が効くのはここ 1 箇所である(**段数分の増幅が起きない**)。

### §2 具体層 —— 根の重複集合と `Γ_K`

* `aroots_map_smul` —— `K` 係数多項式の `K̄` での根は `Γ_K` で**置換される**。
* `exists_algebraMap_eq_aroots_sum` —— ★**同変な関数の根上の和は `K` に入る**。
  ★中間体も跡写像も作らない(`lean-idioms.md` #59/#69 を回避)。
* `exists_algEquiv_smul_eq_of_mem_aroots` / `norm_aeval_eq_of_mem_aroots` ——
  共役どうしはノルムが等しく、`K` 係数多項式の値のノルムも等しい。

### §3 ★★**Sen 型の重み付き評価**(本ファイルの主結果 1)

`exists_norm_sub_algebraMap_le_of_weights` ——
`w : K̄ → K̄` が `Γ_K`-同変で、`x` の共役上で `Σ w(r) = 1`・`‖w(r)‖ ≤ B` なら

  `∃ y ∈ K, ‖x − y‖ ≤ B ε`。

★**重心(`w ≡ 1/n`)はその特別な場合**である
(`exists_norm_sub_algebraMap_le_div_norm_natDegree'` で実際に再導出した)。
★`w(r) = g(r)/f′(r)`(`g ∈ K[X]`、`f = minpoly`)を入れると
`B = ‖g(x)‖/‖f′(x)‖` になり、**`f′(x)` は different**である
(`exists_norm_sub_algebraMap_le_of_polynomial_weights`)。

### §3′ ★★★**different による評価(無条件)**

`aroots_sum_div_derivative_eq_one` —— **Euler の等式** `Σ_r g(r)/f′(r) = 1`
(`g` がモニックで `deg g = n−1` なら成り立つ)。★`Lagrange.coeff_eq_sum` 1 本で出る。
★★**跡写像も中間体も `PowerBasis` も使わない。**

これで §3 の仮定が外れ、次が**無条件**で言える:

| 宣言 | 主張 |
|---|---|
| `exists_norm_sub_algebraMap_le_div_norm_derivative` | `‖x − y‖ ≤ (‖g(x)‖/‖f′(x)‖) ε`(`g` モニック・次数 `n−1`) |
| `exists_norm_sub_algebraMap_le_pow_div_norm_derivative` | `‖x − y‖ ≤ (‖x‖^{n−1}/‖f′(x)‖) ε` |
| ★`exists_norm_sub_algebraMap_le_inv_norm_derivative` | `‖x‖ ≤ 1` なら **`‖x − y‖ ≤ ‖D‖⁻¹ ε`** |

★最後の 1 本が古典的な**「different による Sen の評価」そのもの**である
(`D = f′(x)` は `K(x)/K` の different の生成元)。

### §4 ★★**降下による `AxLemma K C`**(本ファイルの主結果 2)

`AxDescentStep K C` —— 「`[K(x):K] > 1` なら、次数が真に小さく、
`‖x − x′‖ ≤ C ε` で、しかも **`∀σ ‖σx′ − x′‖ ≤ ε`(もとの `ε` を保つ)**
な `x′` が取れる」。

* `axLemma_of_descentStep : AxDescentStep K C → AxLemma K C`
* ★★`axSenTate_of_axDescentStep : AxDescentStep K C → AxSenTate K`

★★**定数が段数によらない**のが要点である。超距離なので段ごとの誤差が
`max` で吸収され、`C` が掛かり算されない。
★★「もとの `ε` を保つ」という条件が**その代償**である ——
これが無いと `ε` が段ごとに `C` 倍され、`C^k` に発散する。
★非空虚性: `descentStep_of_natDegree_tame` —— tame な `x` については
この降下が**実際に取れる**(重心が `K` に落ちるので 1 段で終わる)。

★★**`C = 1` では `AxDescentStep` も `AxLemma` も偽である**(反例は
`AxDescentStep` の docstring。`K = ℚ_p(ζ_p)`, `x = p^{1/p}`)。
☆★本波は最初 `C = 1` 固定で書いてしまい、**書いた後に自分で反例に気づいて**
定数付きに直した。★見立てが外れた点として記録しておく。

## ★★★測って分かったこと —— 「単発の平均化」では Ax の定数は出ない

★§3 の枠組みで `B` を **`x` に依らずに** 取れるか、を測った。答えは **否**である。

`w` が `Γ_K`-同変で `Σ_{r} w(r) = 1` であることは、
`θ := w(x)` が `K(x)` の元で `Tr_{K(x)/K}(θ) = 1` であることと**同じ**である
(同変性から `θ` は `Gal(K̄/K(x))` で固定され、根上の和はちょうど跡)。
そして `Tr(θ) = 1` を満たす `θ` のノルムの最小値は
**different の逆 `‖D_{K(x)/K}‖⁻¹`** であり、これは `x` について**非有界**
(例: `K(p^{1/p^n})` で `v(D) → ∞`)。

⇒ ★★**`∃ B, ∀ x, ∃ w …` は偽である。**
★したがって Ax の定数 `|p|^{-p/(p-1)^2}` は
**「1 回の平均化」では絶対に出ない** —— 塔に沿った**降下**(§4)が要る。
★これは本波で**測って分かった否定的な事実**であり、次の波が
「重み付き重心を頑張れば閉じる」に費やすのを止めるために書いておく。

★★**見立てが外れた点**: 本波の当初の見立ては
「different の評価(`Tr θ = 1` なる `θ` を取る)を一様化すれば `C` が出る」だった。
★**外れた。** 一様化できないことが上の観察で分かった。
★正しい構造は §4 の降下であり、**各段の定数は 1 ではなく `|π|^{-i}`(`i` は分岐の跳び)**、
その総和が収束するというのが Sen の議論の本体である。
本ファイルは**降下の骨格**と**1 段の評価の枠組み**までを用意し、
「1 段の定数を分岐の跳びで評価する」ところを**未着手**として残す。

## ★残った穴(正直に書く)

★★**`AxLemma K C` の項は作れていない。**(仮説なしで `C` を出してはいない。)作れたのは

1. 1 段の評価の一般形(`exists_norm_sub_algebraMap_le_of_weights`)、
2. ★その **different による無条件版**
   (`exists_norm_sub_algebraMap_le_inv_norm_derivative` ほか。Euler の等式込み)、
3. 降下から `AxLemma K C` / `AxSenTate K` を出す骨格
   (`axLemma_of_descentStep` / `axSenTate_of_axDescentStep`)、
4. その降下が tame な `x` については実際に取れること
   (`descentStep_of_natDegree_tame`)、

の 4 つである。★2 の定数 `‖D‖⁻¹` は `x` について**非有界**なので
(上の「測って分かったこと」)、2 だけでは `AxLemma` は出ない。残っているのは
★**「wild(`p ∣ [K(x):K]`)な `x` について、次数を下げつつ
`‖σx′ − x′‖ ≤ ε` を保つ `x′`」の構成**
——これが Sen の補題そのものであり、**数学が足りない**(配管ではない)。
★必要な新ノード:「巡回 `p` 次拡大の跳び `i` と `‖σx − x‖` の関係」
(1 段の定数は `|π_L|^{−i}` になるはず。木の `LowerRamificationGroup.lean` /
`HerbrandFunction.lean` / `HasseArf*.lean` が入力になりうる)。

## 逸脱の記録(CLAUDE.md「逸脱」)

1. `AxDescentStep` は原典 (Ax 1970 / Sen 1972) にこの形では現れない。
   原典は「次数についての帰納法」を地の文で回すが、そこで**定数が段ごとに掛かる**
   のを避けるため、本ファイルは
   **「1 段で `‖σx′ − x′‖ ≤ ε`(もとの `ε`)を保つ」**という条件を足した。
   ★これは強めた仮説であり、その代わり結論の定数が `C` のまま通る。
   ★原典の `C = |p|^{-p/(p-1)^2}` をこの形で得られるかどうかは**未確認**である。
2. `w` の同変性は `K̄` 全体で要求した(共役上だけで十分だが、
   実際に使う重み `g(r)/f′(r)` は全体で同変なので弱める意味がない)。
3. `.src` は `AxLemma.lean` と**同じ項目**(pGC 物理 p.6 Corollary 3.1)を指す。
   本ファイルの中身は `AxSenTate` の入力であって、原典が独立に立てた項目ではない。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped NormedField Valued

/-! ## §1 抽象核

★★以下の 3 本には**分岐・付値・Galois・p 進の語彙が 1 つも出てこない**。 -/

section AbstractCore

/-- **抽象核** —— `Σ w(r)(x − r) = (Σ w(r))·x − Σ w(r)r`。ただの可換環でよい。 -/
theorem multisetSum_weighted_const_sub {A : Type*} [CommRing A]
    (s : Multiset A) (w : A → A) (x : A) :
    (s.map (fun r => w r * (x - r))).sum
      = (s.map w).sum * x - (s.map (fun r => w r * r)).sum := by
  induction s using Multiset.induction with
  | empty => simp
  | cons a t ih =>
      simp only [Multiset.map_cons, Multiset.sum_cons, ih]
      ring

/-- ★★★**抽象核** —— 「重み付き重心」による近似。

`Σ_{r ∈ s} w(r) = 1` かつ `‖w(r)‖ ≤ B` で、各 `r` が `x` から `ε` 以内なら、
重み付き平均 `Σ w(r) r` は `x` から **`B ε` 以内**にある。

★超距離性が効くのはここ 1 箇所である(**項数分の増幅が起きない**)。
★`AxLemma.lean` の `norm_sub_multisetSum_div_le`(重心)は
`w ≡ (card s)⁻¹`・`B = ‖card s‖⁻¹` の場合にあたる。 -/
theorem norm_sub_multisetSum_weighted_le {A : Type*} [NormedField A]
    [IsUltrametricDist A] (s : Multiset A) (w : A → A) {x : A} {B ε : ℝ}
    (hB : 0 ≤ B) (hε : 0 ≤ ε)
    (hsum : (s.map w).sum = 1)
    (hw : ∀ r ∈ s, ‖w r‖ ≤ B)
    (hr : ∀ r ∈ s, ‖x - r‖ ≤ ε) :
    ‖x - (s.map (fun r => w r * r)).sum‖ ≤ B * ε := by
  have key : x - (s.map (fun r => w r * r)).sum = (s.map (fun r => w r * (x - r))).sum := by
    rw [multisetSum_weighted_const_sub, hsum, one_mul]
  rw [key]
  refine norm_multisetSum_le_of_forall_le (by positivity) ?_
  intro a ha
  obtain ⟨r, hr', rfl⟩ := Multiset.mem_map.mp ha
  rw [norm_mul]
  exact mul_le_mul (hw r hr') (hr r hr') (norm_nonneg _) hB

/-- ★★★**抽象核** —— **超距離での降下**。

`deg : M → ℕ` が「複雑さ」で、

* `hbase` —— `deg x = 0` なら `x` は `S` から `ε` 以内、
* `hstep` —— `deg x ≠ 0` なら、`deg` が真に小さく `‖x − x′‖ ≤ ε` の `x′` が
  性質 `P` を保ったまま取れる

とすると、`P` を満たすすべての `x` が `S` から `ε` 以内にある。

★★**要点は「段数によらない」こと**である。超距離なので段ごとの誤差が
`max` で吸収され、`ε` が増幅しない。アルキメデス的なノルムでは成り立たない。
★`M` にも `S` にも構造は要らない(`S` は閉である必要すらない)。
★証明は `deg` についての強帰納法 1 本。 -/
theorem exists_mem_of_descent {M : Type*} [SeminormedAddCommGroup M] [IsUltrametricDist M]
    {S : Set M} (deg : M → ℕ) (P : M → Prop) {ε : ℝ}
    (hbase : ∀ x, P x → deg x = 0 → ∃ y ∈ S, ‖x - y‖ ≤ ε)
    (hstep : ∀ x, P x → deg x ≠ 0 → ∃ x', P x' ∧ deg x' < deg x ∧ ‖x - x'‖ ≤ ε) :
    ∀ x, P x → ∃ y ∈ S, ‖x - y‖ ≤ ε := by
  intro x hx
  generalize hn : deg x = n
  induction n using Nat.strong_induction_on generalizing x with
  | _ n ih =>
    rcases Nat.eq_zero_or_pos n with h0 | hpos
    · exact hbase x hx (hn.trans h0)
    · obtain ⟨x', hx', hlt, hd⟩ := hstep x hx (by omega)
      obtain ⟨y, hyS, hy⟩ := ih (deg x') (by omega) x' hx' rfl
      refine ⟨y, hyS, ?_⟩
      have hmax := IsUltrametricDist.norm_add_le_max (x - x') (x' - y)
      have heq : x - x' + (x' - y) = x - y := by abel
      rw [heq] at hmax
      exact hmax.trans (max_le hd hy)

end AbstractCore

/-! ## §2 具体層 —— 根の重複集合と `Γ_K` -/

variable {p : ℕ} [Fact p.Prime]

/-- ★`K` 係数多項式の `K̄` での根の重複集合は `Γ_K` の作用で**不変**
(根が置換されるだけ)。

`Polynomial.Splits.roots_map_of_injective` の 1 行で出る。 -/
theorem aroots_map_smul (K : PAdicLocalField p) (f : Polynomial K.carrier) (σ : K.absGal) :
    (f.aroots K.closure).map (fun r => σ • r) = f.aroots K.closure := by
  classical
  set i : K.closure →+* K.closure :=
    (σ : K.closure ≃ₐ[K.carrier] K.closure).toRingEquiv.toRingHom with hi
  have hinj : Function.Injective i := σ.injective
  set g : Polynomial K.closure := f.map (algebraMap K.carrier K.closure) with hg
  have hsplits : g.Splits := IsAlgClosed.splits g
  have hcomp : i.comp (algebraMap K.carrier K.closure) = algebraMap K.carrier K.closure :=
    RingHom.ext fun a => σ.commutes a
  have hmap : g.map i = g := by rw [hg, Polynomial.map_map, hcomp]
  have h1 : (g.map i).roots = g.roots.map i := hsplits.roots_map_of_injective hinj
  rw [hmap] at h1
  rw [Polynomial.aroots_def]
  exact h1.symm

def aroots_map_smul.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★**同変な関数の、根の重複集合上の和は `K` に入る**。

★中間体も跡写像も作らない —— `Γ_K` 不変であることを `aroots_map_smul` で示し、
`AxLemma.lean` の `mem_range_of_forall_smul_eq`(`K̄^{Γ_K} = K`)に渡すだけ。
★`lean-idioms.md` #59/#69(中間体の 2 層をまたぐ `rfl`)を完全に回避している。 -/
theorem exists_algebraMap_eq_aroots_sum (K : PAdicLocalField p)
    (f : Polynomial K.carrier) (F : K.closure → K.closure)
    (hequiv : ∀ (σ : K.absGal) (r : K.closure), F (σ • r) = σ • F r) :
    ∃ y : K.carrier, algebraMap K.carrier K.closure y
      = ((f.aroots K.closure).map F).sum := by
  classical
  set s : Multiset K.closure := f.aroots K.closure with hs
  have hinv : ∀ σ : K.absGal, σ • (s.map F).sum = (s.map F).sum := by
    intro σ
    rw [smul_closure_def]
    have h1 : σ ((s.map F).sum) = ((s.map F).map (fun z => σ z)).sum := by
      simpa using map_multiset_sum (σ : K.closure ≃ₐ[K.carrier] K.closure) (s.map F)
    rw [h1, Multiset.map_map]
    have h2 : (fun z => σ z) ∘ F = F ∘ (fun r => σ • r) := by
      funext r
      simp only [Function.comp_apply]
      rw [hequiv σ r, smul_closure_def]
    rw [h2, ← Multiset.map_map F (fun r => σ • r), hs, aroots_map_smul K f σ]
  obtain ⟨y, hy⟩ := mem_range_of_forall_smul_eq K hinv
  exact ⟨y, hy⟩

def exists_algebraMap_eq_aroots_sum.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- `x` の共役は `Γ_K` の軌道である(`IsConjRoot.exists_algEquiv` の言い換え)。 -/
theorem exists_algEquiv_smul_eq_of_mem_aroots (K : PAdicLocalField p) {x r : K.closure}
    (hr : r ∈ (minpoly K.carrier x).aroots K.closure) : ∃ σ : K.absGal, σ • x = r := by
  have hint : IsIntegral K.carrier x :=
    (Algebra.IsAlgebraic.isAlgebraic (R := K.carrier) x).isIntegral
  have hconj : IsConjRoot K.carrier x r := (isConjRoot_iff_mem_minpoly_aroots hint).mpr hr
  obtain ⟨σ, hσ⟩ := hconj.exists_algEquiv
  refine ⟨σ⁻¹, ?_⟩
  rw [← hσ]
  show σ⁻¹ • (σ • r) = r
  rw [inv_smul_smul]

def exists_algEquiv_smul_eq_of_mem_aroots.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★共役における `K` 係数多項式の値のノルムは等しい。

★これが「different `f′(x)` のノルムが共役に依らない」ことの中身である。 -/
theorem norm_aeval_eq_of_mem_aroots (K : PAdicLocalField p) {x r : K.closure}
    (hr : r ∈ (minpoly K.carrier x).aroots K.closure) (h : Polynomial K.carrier) :
    ‖Polynomial.aeval r h‖ = ‖Polynomial.aeval x h‖ := by
  obtain ⟨σ, hσ⟩ := exists_algEquiv_smul_eq_of_mem_aroots K hr
  rw [← hσ, smul_closure_def, Polynomial.aeval_algHom_apply
    (σ : K.closure ≃ₐ[K.carrier] K.closure) x h, ← smul_closure_def, norm_smul_closure]

def norm_aeval_eq_of_mem_aroots.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-! ## §3 主結果 1 —— Sen 型の重み付き評価 -/

/-- ★★★**Sen 型の重み付き評価**。

`w : K̄ → K̄` が `Γ_K`-同変で、`x` の共役の重複集合の上で

* `Σ_r w(r) = 1`、
* `‖w(r)‖ ≤ B`

を満たすなら、`∀σ ‖σx − x‖ ≤ ε` から **`∃ y ∈ K, ‖x − y‖ ≤ B ε`**。

★`w ≡ 1/n` が `AxLemma.lean` の重心にあたる
(`exists_norm_sub_algebraMap_le_div_norm_natDegree'` で再導出した)。
★`w(r) = g(r)/f′(r)` が **different による評価**にあたる
(`exists_norm_sub_algebraMap_le_of_polynomial_weights`)。

★★**注意**: `B` を `x` に依らず取ることは**できない**(モジュール docstring 参照)。 -/
theorem exists_norm_sub_algebraMap_le_of_weights (K : PAdicLocalField p)
    {x : K.closure} {ε B : ℝ} (hε : 0 ≤ ε) (hB : 0 ≤ B)
    (w : K.closure → K.closure)
    (hequiv : ∀ (σ : K.absGal) (r : K.closure), w (σ • r) = σ • w r)
    (hsum : (((minpoly K.carrier x).aroots K.closure).map w).sum = 1)
    (hnorm : ∀ r ∈ (minpoly K.carrier x).aroots K.closure, ‖w r‖ ≤ B)
    (hx : ∀ σ : K.absGal, ‖σ • x - x‖ ≤ ε) :
    ∃ y : K.carrier, ‖x - algebraMap K.carrier K.closure y‖ ≤ B * ε := by
  classical
  have hequiv' : ∀ (σ : K.absGal) (r : K.closure),
      (fun r => w r * r) (σ • r) = σ • (fun r => w r * r) r := by
    intro σ r
    simp only [smul_closure_def]
    rw [map_mul, ← smul_closure_def, hequiv σ r, smul_closure_def]
  obtain ⟨y, hy⟩ :=
    exists_algebraMap_eq_aroots_sum K (minpoly K.carrier x) (fun r => w r * r) hequiv'
  refine ⟨y, ?_⟩
  rw [hy]
  exact norm_sub_multisetSum_weighted_le _ w hB hε hsum hnorm
    (fun r hr => norm_sub_le_of_mem_aroots K hx hr)

def exists_norm_sub_algebraMap_le_of_weights.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★**非空虚性 / 整合性の検査** —— `w ≡ n⁻¹`(重心)を入れると
`AxLemma.lean` の `exists_norm_sub_algebraMap_le_div_norm_natDegree` が再導出できる。

★`exists_norm_sub_algebraMap_le_of_weights` が**空虚でない**ことの証拠であり、
同時に「重心は重み付き重心の特別な場合」であることの確認でもある。 -/
theorem exists_norm_sub_algebraMap_le_div_norm_natDegree'
    (K : PAdicLocalField p) {x : K.closure} {ε : ℝ} (hε : 0 ≤ ε)
    (hx : ∀ σ : K.absGal, ‖σ • x - x‖ ≤ ε) :
    ∃ y : K.carrier, ‖x - algebraMap K.carrier K.closure y‖
      ≤ ‖(((minpoly K.carrier x).natDegree : ℕ) : K.carrier)‖⁻¹ * ε := by
  classical
  set n : ℕ := (minpoly K.carrier x).natDegree with hn
  have hint : IsIntegral K.carrier x :=
    (Algebra.IsAlgebraic.isAlgebraic (R := K.carrier) x).isIntegral
  have hmonic : (minpoly K.carrier x).Monic := minpoly.monic hint
  have hdeg : n ≠ 0 := (minpoly.natDegree_pos hint).ne'
  have hnpos : 0 < ‖((n : ℕ) : K.closure)‖ := norm_natCast_closure_pos K hdeg
  have hne : ((n : ℕ) : K.closure) ≠ 0 := norm_pos_iff.mp hnpos
  have hcard : Multiset.card ((minpoly K.carrier x).aroots K.closure) = n := by
    rw [Polynomial.aroots_def]
    rw [Polynomial.splits_iff_card_roots.mp
      (IsAlgClosed.splits ((minpoly K.carrier x).map (algebraMap K.carrier K.closure))),
      hmonic.natDegree_map]
  refine exists_norm_sub_algebraMap_le_of_weights K hε (by positivity)
    (fun _ => ((n : ℕ) : K.closure)⁻¹) ?_ ?_ ?_ hx
  · intro σ _
    rw [smul_closure_def, map_inv₀, map_natCast]
  · rw [Multiset.map_const', hcard, Multiset.sum_replicate, nsmul_eq_mul]
    exact mul_inv_cancel₀ hne
  · intro r _
    rw [norm_inv, norm_natCast_closure]

def exists_norm_sub_algebraMap_le_div_norm_natDegree'.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★**different による評価**(Sen の補題の「1 段」の形)。

`g ∈ K[X]` を取り、重み `w(r) = g(r)/f′(r)`(`f = minpoly K x`)を使う。
`Σ_r g(r)/f′(r) = 1` が満たされているとき、

  `∃ y ∈ K, ‖x − y‖ ≤ (‖g(x)‖/‖f′(x)‖) · ε`。

★`f′(x)` は `K(x)/K` の **different の生成元**である
(`aeval_derivative_mem_differentIdeal`, `conductor_mul_differentIdeal`)。
★仮定 `Σ_r g(r)/f′(r) = 1` は `g = X^{n−1}` のとき Euler の等式で自動的に成り立つ
(mathlib では `Module.Basis.traceDual_powerBasis_eq` が入口。本ファイルでは
中間体の巻き込みを避けるため**仮定のまま**にしてある)。 -/
theorem exists_norm_sub_algebraMap_le_of_polynomial_weights (K : PAdicLocalField p)
    {x : K.closure} {ε : ℝ} (hε : 0 ≤ ε) (g : Polynomial K.carrier)
    (hsum : (((minpoly K.carrier x).aroots K.closure).map
        (fun r => Polynomial.aeval r g /
          Polynomial.aeval r (Polynomial.derivative (minpoly K.carrier x)))).sum = 1)
    (hx : ∀ σ : K.absGal, ‖σ • x - x‖ ≤ ε) :
    ∃ y : K.carrier, ‖x - algebraMap K.carrier K.closure y‖
      ≤ (‖Polynomial.aeval x g‖ /
          ‖Polynomial.aeval x (Polynomial.derivative (minpoly K.carrier x))‖) * ε := by
  classical
  refine exists_norm_sub_algebraMap_le_of_weights K hε (by positivity) _ ?_ hsum ?_ hx
  · intro σ r
    simp only [smul_closure_def]
    rw [map_div₀, Polynomial.aeval_algHom_apply (σ : K.closure ≃ₐ[K.carrier] K.closure) r g,
      Polynomial.aeval_algHom_apply (σ : K.closure ≃ₐ[K.carrier] K.closure) r
        (Polynomial.derivative (minpoly K.carrier x))]
  · intro r hr
    rw [norm_div, norm_aeval_eq_of_mem_aroots K hr g,
      norm_aeval_eq_of_mem_aroots K hr (Polynomial.derivative (minpoly K.carrier x))]

def exists_norm_sub_algebraMap_le_of_polynomial_weights.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★★**Euler の等式**(`Σ_r g(r)/f′(r) = 1`)—— **`g` がモニックで `deg g = n − 1` なら成り立つ**。

★これは `Lagrange.coeff_eq_sum`(`LinearAlgebra/Lagrange.lean:495`)を
`P = g` に当てただけである:

```
P.coeff (#s − 1) = Σ_{i ∈ s} P.eval (v i) / Π_{j ≠ i} (v i − v j)
```

節点 `v = id` を `f` の根の集合に取れば、`Π_{j ≠ r}(r − j) = f′(r)`
(`Polynomial.eval_multiset_prod_X_sub_C_derivative`)であり、
`g` がモニックで次数 `n−1` なら左辺は `1` である。

★**跡写像も中間体も `PowerBasis` も要らない。** mathlib の
`Module.Basis.traceDual_powerBasis_eq` を経由すると `K(x)` を作る必要が出るが、
`Lagrange.coeff_eq_sum` なら**根の重複集合の上だけ**で済む
(`lean-idioms.md` #59/#69 の回避)。★原典より短い道である。
★根が重複しないこと(`Separable`)は標数 0(`PerfectField.ofCharZero`)から出る。 -/
theorem aroots_sum_div_derivative_eq_one (K : PAdicLocalField p) (x : K.closure)
    (g : Polynomial K.carrier) (hg : g.Monic)
    (hdeg : g.natDegree + 1 = (minpoly K.carrier x).natDegree) :
    (((minpoly K.carrier x).aroots K.closure).map
      (fun r => Polynomial.aeval r g /
        Polynomial.aeval r (Polynomial.derivative (minpoly K.carrier x)))).sum = 1 := by
  classical
  haveI : CharZero K.carrier :=
    charZero_of_injective_algebraMap (algebraMap ℚ_[p] K.carrier).injective
  have hint : IsIntegral K.carrier x :=
    (Algebra.IsAlgebraic.isAlgebraic (R := K.carrier) x).isIntegral
  have hfmonic : (minpoly K.carrier x).Monic := minpoly.monic hint
  have hfsep : (minpoly K.carrier x).Separable := Algebra.IsSeparable.isSeparable K.carrier x
  rw [Polynomial.aroots_def]
  set f : Polynomial K.carrier := minpoly K.carrier x with hf
  set n : ℕ := f.natDegree with hn
  set F : Polynomial K.closure := f.map (algebraMap K.carrier K.closure) with hF
  set G : Polynomial K.closure := g.map (algebraMap K.carrier K.closure) with hG
  have hFmonic : F.Monic := hfmonic.map _
  have hFsplits : F.Splits := IsAlgClosed.splits F
  have hFsep : F.Separable := hfsep.map
  set s : Multiset K.closure := F.roots with hs
  have hnodup : s.Nodup := Polynomial.nodup_roots hFsep
  have hcard : Multiset.card s = n := by
    rw [hs, Polynomial.splits_iff_card_roots.mp hFsplits, hF, hfmonic.natDegree_map]
  set t : Finset K.closure := s.toFinset with ht
  have htval : t.val = s := by rw [ht, Multiset.toFinset_val, Multiset.dedup_eq_self.mpr hnodup]
  have htcard : t.card = n := by rw [Finset.card, htval, hcard]
  have hsum_eq : ∀ h : K.closure → K.closure, (s.map h).sum = ∑ r ∈ t, h r := by
    intro h; rw [← htval]; rfl
  have hdegG : G.degree < (t.card : WithBot ℕ) := by
    rw [htcard]
    refine lt_of_le_of_lt Polynomial.degree_le_natDegree ?_
    rw [hG, hg.natDegree_map]
    exact_mod_cast (by omega : g.natDegree < n)
  have hlag := Lagrange.coeff_eq_sum (v := (id : K.closure → K.closure)) (s := t)
    (Set.injOn_id _) hdegG
  have hGcoeff : G.coeff (t.card - 1) = 1 := by
    rw [htcard]
    have hGd : G.natDegree = n - 1 := by rw [hG, hg.natDegree_map]; omega
    rw [← hGd]
    exact (hg.map (algebraMap K.carrier K.closure)).coeff_natDegree
  have hprod : ∀ r ∈ t, (∏ j ∈ t.erase r, (id r - id j))
      = Polynomial.aeval r (Polynomial.derivative f) := by
    intro r hr
    have hrs : r ∈ s := by rwa [ht, Multiset.mem_toFinset] at hr
    have hFprod : F = (s.map (fun a => Polynomial.X - Polynomial.C a)).prod :=
      hFsplits.eq_prod_roots_of_monic hFmonic
    have h1 := Polynomial.eval_multiset_prod_X_sub_C_derivative (S := s) (r := r) hrs
    rw [← hFprod] at h1
    have h2 : Polynomial.aeval r (Polynomial.derivative f)
        = Polynomial.eval r (Polynomial.derivative F) := by
      rw [hF, Polynomial.derivative_map, Polynomial.eval_map, ← Polynomial.aeval_def]
    rw [h2, h1]
    show ((t.erase r).val.map (fun j => id r - id j)).prod = _
    rw [Finset.erase_val, htval]
    rfl
  calc (s.map (fun r => Polynomial.aeval r g / Polynomial.aeval r (Polynomial.derivative f))).sum
      = ∑ r ∈ t, (Polynomial.aeval r g / Polynomial.aeval r (Polynomial.derivative f)) :=
        hsum_eq _
    _ = ∑ r ∈ t, (G.eval (id r) / ∏ j ∈ t.erase r, (id r - id j)) := by
        refine Finset.sum_congr rfl fun r hr => ?_
        rw [hprod r hr, hG, Polynomial.eval_map, ← Polynomial.aeval_def]
        rfl
    _ = 1 := by rw [← hlag, hGcoeff]

def aroots_sum_div_derivative_eq_one.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★★**different による評価(無条件)** —— `g` がモニックで `deg g = [K(x):K] − 1` なら

  `∃ y ∈ K, ‖x − y‖ ≤ (‖g(x)‖ / ‖f′(x)‖) · ε`   (`f = minpoly K x`)。

★仮定は `g` についてだけで、`x` については**何も仮定していない**。 -/
theorem exists_norm_sub_algebraMap_le_div_norm_derivative (K : PAdicLocalField p)
    {x : K.closure} {ε : ℝ} (hε : 0 ≤ ε) (g : Polynomial K.carrier) (hg : g.Monic)
    (hdeg : g.natDegree + 1 = (minpoly K.carrier x).natDegree)
    (hx : ∀ σ : K.absGal, ‖σ • x - x‖ ≤ ε) :
    ∃ y : K.carrier, ‖x - algebraMap K.carrier K.closure y‖
      ≤ (‖Polynomial.aeval x g‖ /
          ‖Polynomial.aeval x (Polynomial.derivative (minpoly K.carrier x))‖) * ε :=
  exists_norm_sub_algebraMap_le_of_polynomial_weights K hε g
    (aroots_sum_div_derivative_eq_one K x g hg hdeg) hx

def exists_norm_sub_algebraMap_le_div_norm_derivative.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★★**Sen の different 評価(無条件・`g = X^{n−1}`)** ——

  `∃ y ∈ K, ‖x − y‖ ≤ (‖x‖^{n−1} / ‖f′(x)‖) · ε`   (`n = [K(x):K]`)。

★`f′(x)` は `K(x)/K` の different の生成元である。 -/
theorem exists_norm_sub_algebraMap_le_pow_div_norm_derivative (K : PAdicLocalField p)
    {x : K.closure} {ε : ℝ} (hε : 0 ≤ ε)
    (hx : ∀ σ : K.absGal, ‖σ • x - x‖ ≤ ε) :
    ∃ y : K.carrier, ‖x - algebraMap K.carrier K.closure y‖
      ≤ (‖x‖ ^ ((minpoly K.carrier x).natDegree - 1) /
          ‖Polynomial.aeval x (Polynomial.derivative (minpoly K.carrier x))‖) * ε := by
  have hint : IsIntegral K.carrier x :=
    (Algebra.IsAlgebraic.isAlgebraic (R := K.carrier) x).isIntegral
  have hpos : 0 < (minpoly K.carrier x).natDegree := minpoly.natDegree_pos hint
  have h := exists_norm_sub_algebraMap_le_div_norm_derivative K hε
    (Polynomial.X ^ ((minpoly K.carrier x).natDegree - 1)) (Polynomial.monic_X_pow _)
    (by rw [Polynomial.natDegree_X_pow]; omega) hx
  simpa using h

def exists_norm_sub_algebraMap_le_pow_div_norm_derivative.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★**整な `x` については `d(x,K) ≤ ‖D‖⁻¹ ε`**(`D = f′(x)` は different)。

★古典的な「different による Sen の評価」そのものの形である。 -/
theorem exists_norm_sub_algebraMap_le_inv_norm_derivative (K : PAdicLocalField p)
    {x : K.closure} {ε : ℝ} (hε : 0 ≤ ε) (hx1 : ‖x‖ ≤ 1)
    (hx : ∀ σ : K.absGal, ‖σ • x - x‖ ≤ ε) :
    ∃ y : K.carrier, ‖x - algebraMap K.carrier K.closure y‖
      ≤ ‖Polynomial.aeval x (Polynomial.derivative (minpoly K.carrier x))‖⁻¹ * ε := by
  obtain ⟨y, hy⟩ := exists_norm_sub_algebraMap_le_pow_div_norm_derivative K hε hx
  refine ⟨y, hy.trans ?_⟩
  have hpow : ‖x‖ ^ ((minpoly K.carrier x).natDegree - 1) ≤ 1 :=
    pow_le_one₀ (norm_nonneg x) hx1
  have hinv : (0 : ℝ) ≤
      ‖Polynomial.aeval x (Polynomial.derivative (minpoly K.carrier x))‖⁻¹ :=
    inv_nonneg.mpr (norm_nonneg _)
  have hkey : ‖x‖ ^ ((minpoly K.carrier x).natDegree - 1) /
      ‖Polynomial.aeval x (Polynomial.derivative (minpoly K.carrier x))‖
      ≤ ‖Polynomial.aeval x (Polynomial.derivative (minpoly K.carrier x))‖⁻¹ := by
    rw [div_eq_mul_inv]
    nlinarith [pow_nonneg (norm_nonneg x) ((minpoly K.carrier x).natDegree - 1)]
  exact mul_le_mul_of_nonneg_right hkey hε

def exists_norm_sub_algebraMap_le_inv_norm_derivative.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-! ## §4 主結果 2 —— 降下から `AxLemma K C` -/

/-- ★★**降下の 1 段**(本ファイルが原典から切り出した形)。

「`x ∉ K`(= `[K(x):K] > 1`)で `∀σ ‖σx − x‖ ≤ ε` なら、次を満たす `x′` が取れる」:

* `[K(x′):K] < [K(x):K]`(次数が真に下がる)、
* `‖x − x′‖ ≤ C ε`(1 段の誤差)、
* ★★**`∀σ ‖σx′ − x′‖ ≤ ε`** —— **もとの `ε` を保つ**。

★★3 つ目が要点である。これが無いと `ε` が段ごとに `C` 倍され、
`C^k` に発散して**一様な定数が出ない**。★逆にこれさえあれば
超距離性(`exists_mem_of_descent`)で**段数によらず `C ε`** に収まる。

★★**`C = 1` ではこの仮説は偽である。** 反例(手計算。★本ファイルでは形式化していない):
`p ≥ 2`、`K = ℚ_p(ζ_p)`、`x = p^{1/p}`。`σ x = ζ_p^j x` なので
`ε = max_σ ‖σx − x‖ = |p|^{1/(p−1)+1/p}` だが、`K` の値群は `|p|^{ℤ/(p−1)}` で
`1/p ∉ ℤ/(p−1)` だから任意の `a ∈ K` について `‖x − a‖ = max(‖x‖,‖a‖) ≥ |p|^{1/p}`。
よって `d(x,K) = |p|^{−1/(p−1)}·ε > ε`。★したがって `AxLemma K 1` も**偽**である。
★★これは本波で測って分かったことで、`C` を `1` に固定した形で書き始めた
**当初の設計が誤りだった**ことを意味する(見立てが外れた点)。

★原典 (Ax/Sen) の一様定数は `C = |p|^{−p/(p−1)^2}` である。 -/
def AxDescentStep (K : PAdicLocalField p) (C : ℝ) : Prop :=
  ∀ ε : ℝ, 0 ≤ ε → ∀ x : K.closure, (∀ σ : K.absGal, ‖σ • x - x‖ ≤ ε) →
    (minpoly K.carrier x).natDegree ≠ 1 →
      ∃ x' : K.closure,
        (minpoly K.carrier x').natDegree < (minpoly K.carrier x).natDegree ∧
          ‖x - x'‖ ≤ C * ε ∧ ∀ σ : K.absGal, ‖σ • x' - x'‖ ≤ ε

def AxDescentStep.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★★**降下から Ax の補題が出る**(本ファイルの主結果 2)。

`AxDescentStep K C`(1 段の降下)から **`AxLemma K C`** が出る。
★★**定数が段数によらない**のが要点で、これが `exists_mem_of_descent`(超距離)の効き所。
★誤差は `max` で吸収されるので `C` は掛かり算されない。

★`AxLemma K C → AxSenTate K`(`AxLemma.lean` の `axSenTate_of_axLemma`)と繋げば、
**`AxDescentStep K C` ただ 1 つから Ax–Sen–Tate が出る**
(`axSenTate_of_axDescentStep`)。 -/
theorem axLemma_of_descentStep (K : PAdicLocalField p) {C : ℝ} (hC : 0 ≤ C)
    (h : AxDescentStep K C) : AxLemma K C := by
  classical
  intro ε hε x hx
  have hdegpos : ∀ z : K.closure, 0 < (minpoly K.carrier z).natDegree := by
    intro z
    exact minpoly.natDegree_pos (Algebra.IsAlgebraic.isAlgebraic (R := K.carrier) z).isIntegral
  have key : ∀ z : K.closure, (∀ σ : K.absGal, ‖σ • z - z‖ ≤ ε) →
      ∃ y ∈ Set.range (algebraMap K.carrier K.closure), ‖z - y‖ ≤ C * ε := by
    refine exists_mem_of_descent (fun z => (minpoly K.carrier z).natDegree - 1)
      (fun z => ∀ σ : K.absGal, ‖σ • z - z‖ ≤ ε) ?_ ?_
    · intro z _ hz0
      have h1 : (minpoly K.carrier z).natDegree = 1 := by have := hdegpos z; omega
      have hmem : z ∈ (algebraMap K.carrier K.closure).range :=
        minpoly.natDegree_eq_one_iff.mp h1
      obtain ⟨c, hc⟩ := hmem
      exact ⟨z, ⟨c, hc⟩, by simpa using mul_nonneg hC hε⟩
    · intro z hz hz0
      have h1 : (minpoly K.carrier z).natDegree ≠ 1 := by have := hdegpos z; omega
      obtain ⟨z', hlt, hd, hinv⟩ := h ε hε z hz h1
      exact ⟨z', hinv, by have := hdegpos z; have := hdegpos z'; omega, hd⟩
  obtain ⟨y, ⟨c, hc⟩, hy⟩ := key x hx
  exact ⟨c, by rw [hc]; exact hy⟩

def axLemma_of_descentStep.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★★**`AxDescentStep K C` ただ 1 つから Ax–Sen–Tate が出る**。

`AxLemma.lean` の `axSenTate_of_axLemma` と繋いだもの。
★これで pGC §3 の全部(`AxSenTate.lean` の `hodgeTateDim_trivial_zero_eq_one_iff`)が
**「1 段の降下が取れる」ただ 1 点**に落ちた。 -/
theorem axSenTate_of_axDescentStep (K : PAdicLocalField p) {C : ℝ} (hC : 0 ≤ C)
    (h : AxDescentStep K C) : AxSenTate K :=
  axSenTate_of_axLemma hC (axLemma_of_descentStep K hC h)

def axSenTate_of_axDescentStep.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★**非空虚性** —— tame な `x`(`p ∤ [K(x):K]`)については降下が**実際に取れる**。

重心が `K` に落ちるので 1 段で `[K(x′):K] = 1` になり、`x′ ∈ K` は完全不変なので
「もとの `ε` を保つ」条件も自動的に満たされる。
★`AxDescentStep` の仮説が空集合について語っているのではないことの証拠である。
★★**残っているのは wild(`p ∣ [K(x):K]`)な側だけ**である。 -/
theorem descentStep_of_natDegree_tame (K : PAdicLocalField p) {C : ℝ} (hC : 1 ≤ C)
    {x : K.closure} {ε : ℝ} (hε : 0 ≤ ε)
    (hn : ¬ (p ∣ (minpoly K.carrier x).natDegree))
    (hx : ∀ σ : K.absGal, ‖σ • x - x‖ ≤ ε)
    (hne : (minpoly K.carrier x).natDegree ≠ 1) :
    ∃ x' : K.closure,
      (minpoly K.carrier x').natDegree < (minpoly K.carrier x).natDegree ∧
        ‖x - x'‖ ≤ C * ε ∧ ∀ σ : K.absGal, ‖σ • x' - x'‖ ≤ ε := by
  obtain ⟨y, hy⟩ := exists_norm_sub_algebraMap_le_of_natDegree_tame K hε hn hx
  refine ⟨algebraMap K.carrier K.closure y, ?_, hy.trans ?_, ?_⟩
  · have h1 : (minpoly K.carrier (algebraMap K.carrier K.closure y)).natDegree = 1 := by
      rw [minpoly.eq_X_sub_C K.closure y, Polynomial.natDegree_X_sub_C]
    have h2 : 0 < (minpoly K.carrier x).natDegree :=
      minpoly.natDegree_pos (Algebra.IsAlgebraic.isAlgebraic (R := K.carrier) x).isIntegral
    omega
  · nlinarith
  · intro σ
    rw [smul_closure_def, AlgEquiv.commutes, sub_self, norm_zero]
    exact hε

def descentStep_of_natDegree_tame.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-! ## §5 使っている公理の一覧

★どれも `propext` / `Classical.choice` / `Quot.sound` だけで、`sorryAx` は**無い**。 -/

#print axioms multisetSum_weighted_const_sub
#print axioms norm_sub_multisetSum_weighted_le
#print axioms exists_mem_of_descent
#print axioms aroots_map_smul
#print axioms exists_algebraMap_eq_aroots_sum
#print axioms exists_algEquiv_smul_eq_of_mem_aroots
#print axioms norm_aeval_eq_of_mem_aroots
#print axioms exists_norm_sub_algebraMap_le_of_weights
#print axioms exists_norm_sub_algebraMap_le_div_norm_natDegree'
#print axioms exists_norm_sub_algebraMap_le_of_polynomial_weights
#print axioms aroots_sum_div_derivative_eq_one
#print axioms exists_norm_sub_algebraMap_le_div_norm_derivative
#print axioms exists_norm_sub_algebraMap_le_pow_div_norm_derivative
#print axioms exists_norm_sub_algebraMap_le_inv_norm_derivative
#print axioms AxDescentStep
#print axioms axLemma_of_descentStep
#print axioms axSenTate_of_axDescentStep
#print axioms descentStep_of_natDegree_tame

end ABC3.Found.PGC
