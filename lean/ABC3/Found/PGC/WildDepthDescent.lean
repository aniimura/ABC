import ABC3.Meta.Claim
import Mathlib.GroupTheory.Sylow
import Mathlib.GroupTheory.GroupAction.Quotient
import Mathlib.GroupTheory.GroupAction.MultipleTransitivity
import Mathlib.NumberTheory.Padics.PadicVal.Basic
import Mathlib.Analysis.Normed.Group.Ultra
import Mathlib.Analysis.Normed.Field.Basic
import Mathlib.FieldTheory.Galois.Basic
import Mathlib.FieldTheory.Normal.Basic
import Mathlib.FieldTheory.Separable

/-!
# [pGC] wild 深さを 1 段下げる —— ★持ち場が名指しした「中間体」の形は**偽**である

## ★★★冒頭に置く訂正(次の波を止めるため。CLAUDE.md「逸脱」)

★本ファイルに配られた目標は次の形だった:

  `wildDepth K x = k ≥ 2` ⇒ ある中間体 `M`(`K ⊆ M ⊆ K(x)`)が在って `wildDepth M x = k − 1`

★★**これは偽である。**反例を下で形式化した(`not_forall_exists_relIndex_padicValNat_eq`)。

### なぜ偽か —— 3 行の翻訳

`N` を `K(x)/K` の Galois 閉包、`G = Gal(N/K)`、`H = Gal(N/K(x))` とすると

* `[K(x):K] = [G:H]` なので `wildDepth K x = v_p (H.index)`、
* 中間体 `M`(`K ⊆ M ⊆ K(x)`)は `H ≤ H₁ ≤ G` と 1:1(Galois 対応は包含を逆にする)、
* `M ⊆ K(x)` なので `M(x) = K(x)`、したがって
  `wildDepth M x = v_p [K(x):M] = v_p (H.relIndex H₁)`。

⇒ 目標は「`v_p (H.index) = k+1` なら `v_p (H.relIndex H₁) = k` なる `H₁ ≥ H` が在る」。

★★**反例**: `G = A₄ = alternatingGroup (Fin 4)`、`H` = 1 点の固定部分群(位数 3)、`p = 2`。

* `H.index = 4`(`alt4_stabilizer_index`)なので `v₂ = 2 = k+1`、`k = 1`。
* `H` は**極大**である —— `A₄` は 4 点に原始的に作用するから
  (`alternatingGroup.isPreprimitive_of_three_le_card` +
  `MulAction.IsPreprimitive.isCoatom_stabilizer_of_isPreprimitive`)。
* よって `H ≤ H₁` は `H₁ = H`(relIndex 1、`v₂ = 0`)か `H₁ = ⊤`(relIndex 4、`v₂ = 2`)だけ。
  ★**`v₂ = 1` は取れない。**

★★**体としての実現**(これは手計算であり、形式化していない):
`F = ℚ₂(ζ₇)`(不分岐 3 次)の 1 単数群は `ℤ₂[C₃]` と同型なので、
`F^×/(F^×)²` は 2 次元既約 `𝔽₂[C₃]` 加群を含む。対応する `V₄` 拡大 `E/F` は `ℚ₂` 上 Galois で
`Gal(E/ℚ₂) ≅ V₄ ⋊ C₃ = A₄`(`|C₃|` と `|V₄|` は互いに素なので `H²` が消えて分裂する)。
`L = E^{C₃}` は `ℚ₂` の 4 次拡大で**中間体を持たない**。`L = ℚ₂(x)` と書けば
`wildDepth ℚ₂ x = v₂ 4 = 2` だが、深さ 1 の中間体は無い。
★同じことは `(S₄, S₃)` でも起きる(`S₃` は `S₄` で極大、指数 4)。

### ★「`≤ k−1` に弱める」でも駄目

`M = K(x)` 自身を取れば `wildDepth M x = 0 ≤ k−1` なので、
★**中間体の言葉のままでは条件が空虚になる**。`[M:K]` に制御を付けようとすると
上の反例に戻る。⇒ ★**中間体という枠組み自体を捨てるのが正しい。**

## ★★★正しい形 —— `x'` を `K(x)` の外に出す

`AxWildDescent K c`(`Found/PGC/AxTowerDecay.lean:473`)は
★**`x'` が `K(x)` に入ることを要求していない**。この自由度が本質である。正しい降下は:

  `P` を `H` の `p`-Sylow、`Q` を `P ≤ Q ≤ G`・`[Q:P] = p` なるもの(★`Q` は `H` を含まない)、
  `y := (1/p) Σ_{c ∈ Q/P} c • x`(`x` は `P` 不変なので代表元の取り方に依らない)。

このとき

* `y` は `Q` 不変 ⇒ `y ∈ N^Q` かつ `v_p [N^Q : K] = v_p (Q.index) = k − 1`
  ⇒ ★`wildDepth K y ≤ k − 1`、
* `‖y − x‖ ≤ ‖(1/p)‖ · Δ(x)`、
* ★★`‖σ y − y‖ ≤ ‖(1/p)‖ · Δ(x)`(`ε` は `‖1/p‖` 倍しか増えない)。

★指数の勘定はぴったり合う: `P` は `H` の `p`-Sylow なので `v_p(P.index) = v_p(H.index) = k`、
`[Q:P] = p` なので `v_p(Q.index) = k − 1`。★**1 段ちょうど下がる。**
(`exists_pgroup_descent`。★これが `A₄` の反例と両立するのは、
`Q` が `H` を含まない —— つまり `N^Q` が `K(x)` の中間体ではない —— からである。
`A₄` で言えば `H = C₃`、`P = 1`、`Q = C₂ ⊂ V₄` で `[G:Q] = 6`、`v₂ 6 = 1`。★確かに 1 下がる。)

## 在庫の測定(★MCP は 1 度も使っていない。`node tools/leanfile.mjs` のみ)

```
grep -nE "\tIsPGroup\." .cache/mathlib-index.txt
  → Sylow.exists_subgroup_card_pow_succ が在る
    (`p^(n+1) ∣ |G|`・`|H| = p^n` ⇒ `∃ K ⊇ H, |K| = p^(n+1)`)
  ★★これ 1 本で「p 群の中を 1 段ずつ上がる」が済む。自作の鎖は要らなかった。
grep -nE "\tSubgroup\.relindex" .cache/mathlib-index.txt          → ★1 件しか出ない
  ★★綴りが変わっている: 現行 mathlib は `Subgroup.relIndex`(I が大文字)。
  ★`relindex` で 0 件を見て「無い」と書いてはいけない実例。
grep -nE "\tSubgroup\.relIndex" .cache/mathlib-index.txt
  → relIndex_mul_index / relIndex_top_right / relIndex_self / relIndex_eq_one … 30 本
grep -nE "\tMulAction\.index_stabilizer" .cache/mathlib-index.txt
  → index_stabilizer_of_transitive : (stabilizer G x).index = Nat.card X
grep -nE "isPreprimitive" .cache/mathlib-index.txt
  → ★alternatingGroup.isPreprimitive_of_three_le_card(GroupAction/MultipleTransitivity.lean:660)
  → ★MulAction.IsPreprimitive.isCoatom_stabilizer_of_isPreprimitive(Primitive.lean:271)
  ★★この 2 本で `A₄` の反例が証明 6 行で済んだ(自前で `A₄` の部分群束を数えずに済む)。
grep -n "norm_sum_le_of_forall_le" .cache/mathlib-index.txt       → ★0 件
  ★しかし `#check @IsUltrametricDist.norm_sum_le_of_forall_le_of_nonneg` は**在る**。
  ★★「索引に無い ⇒ 不在」の反例がこれで 7 例目である(`AxEpsilonDecay` の 6 例目に続く)。
grep -n "IsUltrametricDist.norm_sub_le_max"                       → ★本当に無い
  (`norm_add_le_max` は在る。`to_additive` が片方しか作っていない)
  ⇒ `sub_eq_add_neg` + `norm_neg` + `norm_add_le_max` で書く。
grep -cE "\tABC3\.[A-Za-z.]*\.<名前>(\.|\t)" .cache/decl-index.txt → 本ファイルの全宣言名で 0
```

## 何が言えたか

### §1 反例(★分岐・付値・Galois・`p` 進の語彙が 1 語も出てこない)

| 宣言 | 内容 |
|---|---|
| `relIndex_eq_one_or_index_of_isCoatom` | 極大部分群の上は `1` か `index` しかない |
| `alt4_stabilizer_isCoatom` / `alt4_stabilizer_index` | `A₄` の点固定群は極大・指数 4 |
| ★★`not_forall_exists_relIndex_padicValNat_eq` | ★**持ち場の主張の否定** |

### §2 正しい降下の核(★同じく群論だけ)

| 宣言 | 内容 |
|---|---|
| `exists_le_card_eq_prime_mul` | `p` 群 `P` の指数が `p` で割れる ⇒ `P ≤ Q` かつ `|Q| = p·|P|` |
| `index_eq_prime_mul_index` / `relIndex_eq_prime` | `P.index = p · Q.index`、`P.relIndex Q = p` |
| `padicValNat_index_eq_succ` | ★`v_p(P.index) = v_p(Q.index) + 1` |
| `exists_isPGroup_le_padicValNat_index_eq` | `H` の `p`-Sylow は `v_p(index)` を変えない |
| ★★`exists_pgroup_descent` | ★`v_p(H.index) = k+1` ⇒ `P ≤ H`・`P ≤ Q`・`[Q:P] = p`・`v_p(Q.index) = k` |

### §3 平均化の核(★超距離だけ。付値も Galois も出てこない)

| 宣言 | 内容 |
|---|---|
| `norm_sum_smul_sub_nsmul_le` | `‖Σ_{g∈s} g•x − |s|•x‖ ≤ Δ` |
| `norm_smul_sum_smul_sub_le` | ★`‖τ•(Σ g•x) − Σ g•x‖ ≤ Δ`(★**和を取っても `ε` は増えない**) |
| `cosetSmul` / `cosetSmul_mk` | `x` が `P` 不変なら `Γ ⧸ P → M`, `mk g ↦ g • x` が定まる |
| ★`smul_sum_quotient_eq` | ★**剰余類上の和は `Γ` 不変** |
| `norm_sum_quotient_sub_nsmul_le` | `‖Σ_c f c − [Γ:P]•x‖ ≤ Δ` |
| `norm_map_sum_quotient_sub_le` | 外から来る加法写像 `φ` でも `‖φ y − y‖ ≤ Δ` |

### §4 抽象核と平均化を貼り合わせた「1 段の降下」

| 宣言 | 内容 |
|---|---|
| ★★★`exists_smul_invariant_of_padicValNat_index_succ` | ★**降下 1 段が抽象群で閉じた** |

### §5 / §6 ★★抽象 Galois 版 —— `wildDepth` の言葉のまま

| 宣言 | 内容 |
|---|---|
| ★★`index_stabilizer_eq_natDegree_minpoly` | ★`[Gal(E/F) : Stab(y)] = deg minpoly_F y` |
| `padicValNat_le_of_dvd` | `a ∣ b` ⇒ `v_p a ≤ v_p b` |
| ★★★`exists_natDegree_minpoly_descent` | ★`v_p(deg minpoly x) = k+1` ⇒ `v_p(deg minpoly y) ≤ k` |
| ★★★`exists_natDegree_minpoly_descent_div` | ★★**`AxWildDescent` と同じ形**(`p` で割った版) |

★★★**`k ≥ 2` と `k = 1` を分ける必要は無かった。**§5 の 1 本が `k ≥ 1` を全部扱う。
★持ち場が「`k ≥ 2` は中間体を作る別の話」と想定していたのは、
★**中間体で降りようとしたためである**(それが偽なのは §1)。

★★**`IntermediateField` を 1 度も作っていない** —— 部分体ではなく
`MulAction.stabilizer` で `wildDepth` を測ったので、`lean-idioms.md` #59 の
「中間体 2 層の `rfl` が kernel を止める」に**一度も触らずに済んだ**。

すなわち: 有限群 `G` が超距離加法群 `M` に作用し、`x` が `H` で固定され `∀g ‖g•x − x‖ ≤ Δ`、
`v_p(H.index) = k+1` なら、`Q ≤ G` と `y ∈ M` が在って

* `v_p(Q.index) = k`、
* `∀ q ∈ Q, q • y = y`、
* `‖y − p • x‖ ≤ Δ`(★`p•x` は `p` 倍。体では両辺を `p` で割る)、
* `∀ g : G, ‖g • y − y‖ ≤ Δ`。

★★**`Δ` が増えていない**ことに注意。伸びるのは「`p` で割る」ところだけである。

## ★残った穴(正直に書く。★ここから先は配管であって数学ではない)

★★**`AxWildDescent K (axDecay p)` の項は作れていない。**本ファイルは抽象 Galois の段で止めた。
☆★**追記(同じ波の中で片付いた)**: 下の 1・2 は
`Found/PGC/WildDepthFieldDescent.lean` で**閉じた**
(`axWildDescent_normInv` / `axWildDescent_prime`、`sorry` 0)。
残っているのは **3 だけ**、すなわち定数を `p` から `p^{(1/(p−1))p^{1−k}}` に絞ることである。

体の層(`PAdicLocalField`)に落とすのに要ったのは次の 3 点で、どれも数学ではなく配管だった。

1. `K.closure` の中に `x` を含む有限次 Galois 中間体 `M` を取り、`G := M ≃ₐ[K.carrier] M` を使う。
   ★★**測った**(`node tools/leanfile.mjs`、2026-09-08):
   `lean-idioms.md` #153 の `FiniteGaloisIntermediateField.adjoin K.carrier ({x} : Set K.closure)`
   は `haveI := isGalois_closure K`(`Found/PGC/SubgroupCorrespondenceConstruction.lean:50`)
   を先に置けば `FiniteDimensional` / `Normal` / `Finite (M ≃ₐ[K.carrier] M)` が
   ★**3 つとも `inferInstance` で出る**。さらに
   `NormedField ↥M` / `IsUltrametricDist ↥M` / `DistribMulAction (M ≃ₐ[K.carrier] M) ↥M`
   も ★**3 つとも `inferInstance`**。⇒ §5・§6 の仮定は**全部揃う**。
   ★★**前の波の「中間体 2 層で越えられない」という判断はここでは当たらない。**
2. `K.absGal` の作用を `G` の作用に落とす(`AlgEquiv.restrictNormalHom` の全射性 +
   `AlgEquiv.restrictNormal_commutes`)、および `IntermediateField.minpoly_eq` による
   `minpoly K.carrier (⟨x, _⟩ : ↥M) = minpoly K.carrier x` の移送。
   ☆★**閉じた**(`WildDepthFieldDescent.lean`)。`‖(a : ↥M)‖ = ‖(a : K.closure)‖` と
   `((a − b : ↥M) : K.closure) = (a : K.closure) − (b : K.closure)` が**どちらも `rfl`** で、
   ★往復は 1 発だった。
3. 得られる定数は ★`‖(p : K.closure)‖⁻¹ = p` であって
   `axDecay p k = p^{(1/(p−1))p^{1−k}}` ではない。
   ★`p^{1/(p−1)}` にするには分岐(異なるイデアル)が要る ——
   それは `Found/PGC/CyclicJumpNorm.lean` の `k = 1` の段(別の agent の持ち場)と同じ話である。
   ★★**`c k ≡ p` では `∏_{k≤n} c k = p^n` が発散するので `AxLemma` は出ない。**
   本ファイルが与えるのは「1 段が閉じる」ことだけで、定数の収束は与えない。
   ☆★ただし `AxTowerDecay.axWildDescent_pow`(`c k = p^k`)よりは真に良い。

## 逸脱の記録(CLAUDE.md「逸脱」)

1. ★★持ち場が名指しした「中間体で 1 段下がる」を**偽と判定して捨てた**。
   代わりに「`x'` を `K(x)` の外に出す」形に置き換えた。理由は冒頭。
2. §2〜§4 は原典 (Ax 1970) に対応する文が無い。原典は定数の勘定を地の文で畳んでおり、
   ここで切り出した群論の段は原典が明示していない。★`.src` は
   `AxLemma.lean` / `AxTowerDecay.lean` / `AxEpsilonDecay.lean` と同じ項目
   (pGC 物理 p.6 Corollary 3.1)を指す。原典が独立に立てた項目ではない。
3. §3 の抽象核は `M` に超距離な半ノルム加法群しか要求せず、作用は `DistribMulAction` でよい
   (★等長性も忠実性も要らない)。原典より弱い設定である。
-/

namespace ABC3.Found.PGC

/-! ## §1 ★持ち場の主張の反例(純群論) -/

section Counterexample

variable {G : Type*} [Group G]

/-- 極大部分群 `H` の上には `H` と `⊤` しか無い。 -/
theorem relIndex_eq_one_or_index_of_isCoatom {H : Subgroup G} (hH : IsCoatom H)
    {H₁ : Subgroup G} (h : H ≤ H₁) :
    H.relIndex H₁ = 1 ∨ H.relIndex H₁ = H.index := by
  rcases eq_or_lt_of_le h with rfl | hlt
  · exact Or.inl (Subgroup.relIndex_self _)
  · rw [hH.2 _ hlt, Subgroup.relIndex_top_right]
    exact Or.inr rfl

def relIndex_eq_one_or_index_of_isCoatom.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- `A₄` の 1 点固定部分群(位数 3)は極大。

★4 点への作用が原始的であることから出る。`A₄` の部分群束は数えていない。 -/
theorem alt4_stabilizer_isCoatom :
    IsCoatom (MulAction.stabilizer (alternatingGroup (Fin 4)) (0 : Fin 4)) := by
  haveI : MulAction.IsPreprimitive (alternatingGroup (Fin 4)) (Fin 4) :=
    alternatingGroup.isPreprimitive_of_three_le_card (Fin 4) (by simp)
  exact MulAction.IsPreprimitive.isCoatom_stabilizer_of_isPreprimitive _ _

def alt4_stabilizer_isCoatom.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- `[A₄ : C₃] = 4`。 -/
theorem alt4_stabilizer_index :
    (MulAction.stabilizer (alternatingGroup (Fin 4)) (0 : Fin 4)).index = 4 := by
  haveI : MulAction.IsPreprimitive (alternatingGroup (Fin 4)) (Fin 4) :=
    alternatingGroup.isPreprimitive_of_three_le_card (Fin 4) (by simp)
  rw [MulAction.index_stabilizer_of_transitive]
  simp

def alt4_stabilizer_index.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- `v₂ [A₄ : C₃] = 2` —— ★深さ `2` の状況が実在する。 -/
theorem alt4_padicValNat_index :
    padicValNat 2 (MulAction.stabilizer (alternatingGroup (Fin 4)) (0 : Fin 4)).index = 1 + 1 := by
  haveI : Fact (Nat.Prime 2) := ⟨Nat.prime_two⟩
  rw [alt4_stabilizer_index, show (4 : ℕ) = 2 ^ 2 by norm_num, padicValNat.prime_pow]

def alt4_padicValNat_index.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★`A₄ ⊃ C₃` には「深さがちょうど 1 下がる中間部分群」が無い。 -/
theorem alt4_no_intermediate
    (H₁ : Subgroup (alternatingGroup (Fin 4)))
    (h : MulAction.stabilizer (alternatingGroup (Fin 4)) (0 : Fin 4) ≤ H₁) :
    padicValNat 2 ((MulAction.stabilizer (alternatingGroup (Fin 4)) (0 : Fin 4)).relIndex H₁)
      ≠ 1 := by
  haveI : Fact (Nat.Prime 2) := ⟨Nat.prime_two⟩
  rcases relIndex_eq_one_or_index_of_isCoatom alt4_stabilizer_isCoatom h with h1 | h1
  · rw [h1, padicValNat_one_right]; omega
  · rw [h1, alt4_padicValNat_index]; omega

def alt4_no_intermediate.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★★**持ち場が名指しした主張は偽である。**

Galois 対応で言い換えると、目標は
「`v_p (H.index) = k+1` なら `H ≤ H₁` で `v_p (H.relIndex H₁) = k` なるものが在る」だった。
★`G = A₄`、`H` = 1 点固定部分群、`p = 2`、`k = 1` が反例。 -/
theorem not_forall_exists_relIndex_padicValNat_eq :
    ¬ ∀ (G : Type) [Group G] [Finite G] (p : ℕ) (_ : p.Prime) (H : Subgroup G) (k : ℕ),
        padicValNat p H.index = k + 1 →
        ∃ H₁ : Subgroup G, H ≤ H₁ ∧ padicValNat p (H.relIndex H₁) = k := by
  intro h
  obtain ⟨H₁, hle, hv⟩ :=
    h (alternatingGroup (Fin 4)) 2 Nat.prime_two
      (MulAction.stabilizer (alternatingGroup (Fin 4)) (0 : Fin 4)) 1 alt4_padicValNat_index
  exact alt4_no_intermediate H₁ hle hv

def not_forall_exists_relIndex_padicValNat_eq.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

end Counterexample

/-! ## §2 ★正しい降下の核(純群論) -/

section PGroupDescent

variable {G : Type*} [Group G] [Finite G] {p : ℕ} [hp : Fact p.Prime]

omit hp in
/-- 有限群では指数は `0` でない。 -/
theorem index_ne_zero_of_finite (H : Subgroup G) : H.index ≠ 0 :=
  Subgroup.finiteIndex_iff.mp inferInstance

def index_ne_zero_of_finite.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★`p` 部分群 `P` の指数が `p` で割れるなら、`P` を含む `Q` で `|Q| = p·|P|` なるものが在る。

★mathlib の `Sylow.exists_subgroup_card_pow_succ` 1 本で済む。 -/
theorem exists_le_card_eq_prime_mul {P : Subgroup G} (hP : IsPGroup p P)
    (hdvd : p ∣ P.index) : ∃ Q : Subgroup G, P ≤ Q ∧ Nat.card Q = p * Nat.card P := by
  obtain ⟨n, hn⟩ := IsPGroup.iff_card.mp hP
  obtain ⟨m, hm⟩ := hdvd
  have hcard : Nat.card G = p ^ (n + 1) * m := by
    rw [← Subgroup.card_mul_index P, hn, hm]; ring
  obtain ⟨K, hK, hPK⟩ := Sylow.exists_subgroup_card_pow_succ (G := G) (p := p) (n := n)
      (by rw [hcard]; exact Dvd.intro m rfl) hn
  exact ⟨K, hPK, by rw [hK, hn]; ring⟩

def exists_le_card_eq_prime_mul.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

omit hp in
/-- `|Q| = p·|P|` なら `[G:P] = p·[G:Q]`。 -/
theorem index_eq_prime_mul_index {P Q : Subgroup G} (h : Nat.card Q = p * Nat.card P) :
    P.index = p * Q.index := by
  have h1 : Nat.card P * P.index = Nat.card G := Subgroup.card_mul_index P
  have h2 : Nat.card Q * Q.index = Nat.card G := Subgroup.card_mul_index Q
  have hpos : 0 < Nat.card P := Nat.card_pos
  refine Nat.eq_of_mul_eq_mul_left hpos ?_
  calc Nat.card P * P.index = Nat.card Q * Q.index := by rw [h1, h2]
    _ = Nat.card P * (p * Q.index) := by rw [h]; ring

def index_eq_prime_mul_index.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

omit hp in
/-- `[Q:P] = p`。 -/
theorem relIndex_eq_prime {P Q : Subgroup G} (hPQ : P ≤ Q) (h : Nat.card Q = p * Nat.card P) :
    P.relIndex Q = p := by
  have hmul := Subgroup.relIndex_mul_index hPQ
  rw [index_eq_prime_mul_index h] at hmul
  refine Nat.eq_of_mul_eq_mul_right (Nat.pos_of_ne_zero (index_ne_zero_of_finite Q)) ?_
  rw [hmul]

def relIndex_eq_prime.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★**深さがちょうど 1 下がる**: `v_p [G:P] = v_p [G:Q] + 1`。 -/
theorem padicValNat_index_eq_succ {P Q : Subgroup G} (h : Nat.card Q = p * Nat.card P) :
    padicValNat p P.index = padicValNat p Q.index + 1 := by
  rw [index_eq_prime_mul_index h, padicValNat.mul hp.out.pos.ne' (index_ne_zero_of_finite Q),
    padicValNat.self hp.out.one_lt]
  omega

def padicValNat_index_eq_succ.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★`H` の `p`-Sylow へ落としても `v_p (index)` は変わらない。 -/
theorem exists_isPGroup_le_padicValNat_index_eq (H : Subgroup G) :
    ∃ P : Subgroup G, P ≤ H ∧ IsPGroup p P ∧ padicValNat p P.index = padicValNat p H.index := by
  obtain ⟨S⟩ := (Sylow.nonempty : Nonempty (Sylow p H))
  refine ⟨(S : Subgroup H).map H.subtype, Subgroup.map_subtype_le _, S.isPGroup'.map _, ?_⟩
  rw [Subgroup.index_map_subtype,
    padicValNat.mul (index_ne_zero_of_finite _) (index_ne_zero_of_finite H),
    padicValNat.eq_zero_of_not_dvd S.not_dvd_index]
  omega

def exists_isPGroup_le_padicValNat_index_eq.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★★**正しい 1 段の降下(群論の段)**。

`v_p (H.index) = k+1` なら、`P ≤ H`(`H` の `p`-Sylow)と `P ≤ Q` が在って

* `P` は `p` 群、`[Q:P] = p`、
* `v_p (P.index) = k+1`、★`v_p (Q.index) = k`。

★★`Q` は `H` を含まない —— そこが `A₄` の反例(§1)との違いである。 -/
theorem exists_pgroup_descent (H : Subgroup G) (k : ℕ)
    (hk : padicValNat p H.index = k + 1) :
    ∃ P Q : Subgroup G, P ≤ H ∧ IsPGroup p P ∧ P ≤ Q ∧ P.relIndex Q = p ∧
      padicValNat p P.index = k + 1 ∧ padicValNat p Q.index = k := by
  obtain ⟨P, hPH, hPp, hPv⟩ := exists_isPGroup_le_padicValNat_index_eq (p := p) H
  rw [hk] at hPv
  have hdvd : p ∣ P.index := by
    by_contra hcon
    rw [padicValNat.eq_zero_of_not_dvd hcon] at hPv
    omega
  obtain ⟨Q, hPQ, hQ⟩ := exists_le_card_eq_prime_mul hPp hdvd
  refine ⟨P, Q, hPH, hPp, hPQ, relIndex_eq_prime hPQ hQ, hPv, ?_⟩
  have := padicValNat_index_eq_succ (p := p) hQ
  omega

def exists_pgroup_descent.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

end PGroupDescent

/-! ## §3 ★平均化の核(超距離。★分岐・付値・Galois・`p` 進の語彙が 1 語も出てこない) -/

section Averaging

variable {Γ : Type*} [Group Γ] {M : Type*} [SeminormedAddCommGroup M] [IsUltrametricDist M]
  [DistribMulAction Γ M]

/-- `‖Σ_{g∈s} g•x − |s|•x‖ ≤ Δ`。 -/
theorem norm_sum_smul_sub_nsmul_le {s : Finset Γ} {x : M} {Δ : ℝ} (hΔ : 0 ≤ Δ)
    (hx : ∀ g : Γ, ‖g • x - x‖ ≤ Δ) :
    ‖(∑ g ∈ s, g • x) - s.card • x‖ ≤ Δ := by
  have h : (∑ g ∈ s, g • x) - s.card • x = ∑ g ∈ s, (g • x - x) := by
    rw [Finset.sum_sub_distrib, Finset.sum_const]
  rw [h]
  exact IsUltrametricDist.norm_sum_le_of_forall_le_of_nonneg hΔ (fun g _ => hx g)

def norm_sum_smul_sub_nsmul_le.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★**和を取っても「動きの大きさ」は増えない**。 -/
theorem norm_smul_sum_smul_sub_le {s : Finset Γ} {x : M} {Δ : ℝ} (hΔ : 0 ≤ Δ)
    (hx : ∀ g : Γ, ‖g • x - x‖ ≤ Δ) (τ : Γ) :
    ‖τ • (∑ g ∈ s, g • x) - (∑ g ∈ s, g • x)‖ ≤ Δ := by
  rw [Finset.smul_sum, ← Finset.sum_sub_distrib]
  refine IsUltrametricDist.norm_sum_le_of_forall_le_of_nonneg hΔ (fun g _ => ?_)
  have h : τ • g • x - g • x = ((τ * g) • x - x) - (g • x - x) := by rw [mul_smul]; abel
  rw [h, sub_eq_add_neg]
  refine le_trans (IsUltrametricDist.norm_add_le_max _ _) ?_
  rw [norm_neg]
  exact max_le (hx _) (hx g)

def norm_smul_sum_smul_sub_le.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

variable (P : Subgroup Γ)

/-- `x` が `P` で固定されるとき、剰余類 `Γ ⧸ P` 上の写像 `mk g ↦ g • x`。 -/
def cosetSmul (x : M) (hx : ∀ h ∈ P, h • x = x) : Γ ⧸ P → M := fun c =>
  Quotient.liftOn' c (fun g => g • x) (by
    intro a b hab
    have hmem : a⁻¹ * b ∈ P := QuotientGroup.leftRel_apply.mp hab
    calc a • x = a • ((a⁻¹ * b) • x) := by rw [hx _ hmem]
      _ = b • x := by rw [← mul_smul]; congr 1; group)

def cosetSmul.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

omit [IsUltrametricDist M] in
@[simp] theorem cosetSmul_mk {x : M} (hx : ∀ h ∈ P, h • x = x) (g : Γ) :
    cosetSmul P x hx (QuotientGroup.mk g) = g • x := rfl

def cosetSmul_mk.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

variable [Fintype (Γ ⧸ P)]

omit [IsUltrametricDist M] in
/-- ★★**剰余類上の和は `Γ` 全体で不変**。 -/
theorem smul_sum_quotient_eq {x : M} {f : Γ ⧸ P → M}
    (hf : ∀ g : Γ, f (QuotientGroup.mk g) = g • x) (q : Γ) :
    q • (∑ c : Γ ⧸ P, f c) = ∑ c : Γ ⧸ P, f c := by
  rw [Finset.smul_sum]
  have step : ∀ c : Γ ⧸ P, q • f c = f (q • c) := by
    intro c
    induction c using QuotientGroup.induction_on with
    | _ a =>
      rw [hf a, show q • (QuotientGroup.mk a : Γ ⧸ P) = QuotientGroup.mk (q * a) from rfl,
        hf (q * a), mul_smul]
  simp_rw [step]
  exact Equiv.sum_comp (MulAction.toPerm q) f

def smul_sum_quotient_eq.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- `‖Σ_c f c − [Γ:P]•x‖ ≤ Δ`。 -/
theorem norm_sum_quotient_sub_nsmul_le {x : M} {f : Γ ⧸ P → M} {Δ : ℝ} (hΔ : 0 ≤ Δ)
    (hf : ∀ g : Γ, f (QuotientGroup.mk g) = g • x)
    (hx : ∀ g : Γ, ‖g • x - x‖ ≤ Δ) :
    ‖(∑ c : Γ ⧸ P, f c) - P.index • x‖ ≤ Δ := by
  have hcard : P.index = Finset.univ.card (α := Γ ⧸ P) := by
    rw [Finset.card_univ, Subgroup.index, Nat.card_eq_fintype_card]
  have h : (∑ c : Γ ⧸ P, f c) - P.index • x = ∑ c : Γ ⧸ P, (f c - x) := by
    rw [Finset.sum_sub_distrib, Finset.sum_const, hcard]
  rw [h]
  refine IsUltrametricDist.norm_sum_le_of_forall_le_of_nonneg hΔ (fun c _ => ?_)
  induction c using QuotientGroup.induction_on with
  | _ a => rw [hf a]; exact hx a

def norm_sum_quotient_sub_nsmul_le.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★外から来る加法写像 `φ` に対しても `‖φ y − y‖ ≤ Δ`。

★`φ` は群作用でなくてよい。体の層では `φ = (σ • ·)`(`σ` は大きい群の元)を代入する。 -/
theorem norm_map_sum_quotient_sub_le (φ : M →+ M) {x : M} {f : Γ ⧸ P → M} {Δ : ℝ} (hΔ : 0 ≤ Δ)
    (hf : ∀ g : Γ, f (QuotientGroup.mk g) = g • x)
    (hφ : ∀ g : Γ, ‖φ (g • x) - g • x‖ ≤ Δ) :
    ‖φ (∑ c : Γ ⧸ P, f c) - (∑ c : Γ ⧸ P, f c)‖ ≤ Δ := by
  rw [map_sum, ← Finset.sum_sub_distrib]
  refine IsUltrametricDist.norm_sum_le_of_forall_le_of_nonneg hΔ (fun c _ => ?_)
  induction c using QuotientGroup.induction_on with
  | _ a => rw [hf a]; exact hφ a

def norm_map_sum_quotient_sub_le.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

end Averaging

/-! ## §4 ★★貼り合わせ —— 抽象群における「1 段の降下」 -/

/-- ★★★**降下の 1 段が抽象群で閉じた**(本ファイルの主結果)。

有限群 `G` が超距離な半ノルム加法群 `M` に作用し、`x` は `H ≤ G` で固定され、
`∀ g, ‖g•x − x‖ ≤ Δ` とする。`v_p (H.index) = k+1` なら、`Q ≤ G` と `y ∈ M` が在って

* `v_p (Q.index) = k`(★**深さがちょうど 1 下がる**)、
* `∀ q ∈ Q, q • y = y`、
* `‖y − p • x‖ ≤ Δ`、
* ★`∀ g : G, ‖g • y − y‖ ≤ Δ`(★**`Δ` は増えない**)。

★★体の層では `y/p` を取るので、損失は `‖(p : K)‖⁻¹` ちょうどである。
★`Q` は `H` を含まないので、これは「中間体」ではない。§1 の反例と両立する。 -/
theorem exists_smul_invariant_of_padicValNat_index_succ
    {G : Type*} [Group G] [Finite G] {p : ℕ} [Fact p.Prime]
    {M : Type*} [SeminormedAddCommGroup M] [IsUltrametricDist M] [DistribMulAction G M]
    (H : Subgroup G) {x : M} (hxH : ∀ h ∈ H, h • x = x)
    {Δ : ℝ} (hΔ : 0 ≤ Δ) (hx : ∀ g : G, ‖g • x - x‖ ≤ Δ)
    (k : ℕ) (hk : padicValNat p H.index = k + 1) :
    ∃ (Q : Subgroup G) (y : M), padicValNat p Q.index = k ∧
      (∀ q ∈ Q, q • y = y) ∧ ‖y - p • x‖ ≤ Δ ∧ ∀ g : G, ‖g • y - y‖ ≤ Δ := by
  classical
  obtain ⟨P, Q, hPH, _hPp, hPQ, hrel, _hPv, hQv⟩ := exists_pgroup_descent (p := p) H k hk
  have hxP' : ∀ h ∈ P.subgroupOf Q, (h : ↥Q) • x = x := fun h hh => hxH (h : G) (hPH hh)
  haveI : Fintype (↥Q ⧸ P.subgroupOf Q) := Fintype.ofFinite _
  set f : ↥Q ⧸ P.subgroupOf Q → M := cosetSmul (P.subgroupOf Q) x hxP' with hfdef
  have hf : ∀ g : ↥Q, f (QuotientGroup.mk g) = g • x := fun _ => rfl
  have hxQ : ∀ g : ↥Q, ‖g • x - x‖ ≤ Δ := fun g => hx (g : G)
  refine ⟨Q, ∑ c : ↥Q ⧸ P.subgroupOf Q, f c, hQv, ?_, ?_, ?_⟩
  · intro q hq
    exact smul_sum_quotient_eq (P.subgroupOf Q) hf ⟨q, hq⟩
  · have h := norm_sum_quotient_sub_nsmul_le (P.subgroupOf Q) hΔ hf hxQ
    rwa [show (P.subgroupOf Q).index = p from hrel] at h
  · intro g
    refine norm_map_sum_quotient_sub_le (P.subgroupOf Q)
      (DistribSMul.toAddMonoidHom M g) hΔ hf (fun q => ?_)
    have h : g • ((q : G) • x) - (q : G) • x
        = ((g * (q : G)) • x - x) - ((q : G) • x - x) := by rw [mul_smul]; abel
    show ‖g • ((q : ↥Q) • x) - (q : ↥Q) • x‖ ≤ Δ
    rw [show ((q : ↥Q) • x) = (q : G) • x from rfl, h, sub_eq_add_neg]
    refine le_trans (IsUltrametricDist.norm_add_le_max _ _) ?_
    rw [norm_neg]
    exact max_le (hx _) (hx _)

def exists_smul_invariant_of_padicValNat_index_succ.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-! ## §5 ★★★抽象 Galois 版 —— `wildDepth` の言葉のままで 1 段下げる

★ここでも**分岐・付値・`p` 進の語彙は 1 つも出てこない**(Galois は出る)。
★★`IntermediateField` を 1 度も作らないので、`lean-idioms.md` #59 の
「中間体 2 層の `rfl` が kernel を止める」に**一度も触らない** —— 部分体ではなく
`MulAction.stabilizer` を使うのが鍵である。 -/

/-- ★★**安定化群の指数は最小多項式の次数**(`[Gal(E/F) : Stab(y)] = deg minpoly_F y`)。

★軌道 = `minpoly` の根の集合(`Normal.minpoly_eq_iff_mem_orbit`)+
根の個数 = 次数(`Polynomial.card_rootSet_eq_natDegree`)。
★★これが `wildDepth` を群論に翻訳する橋である:
`[K(x):K] = deg minpoly` なので `wildDepth = v_p (Stab(x).index)`。 -/
theorem index_stabilizer_eq_natDegree_minpoly
    {F : Type*} [Field F] {E : Type*} [Field E] [Algebra F E]
    [FiniteDimensional F E] [hN : Normal F E] [Algebra.IsSeparable F E] (y : E) :
    (MulAction.stabilizer (E ≃ₐ[F] E) y).index = (minpoly F y).natDegree := by
  classical
  have hint : IsIntegral F y := Algebra.IsIntegral.isIntegral y
  have hmonic : (minpoly F y).Monic := minpoly.monic hint
  have hirr : Irreducible (minpoly F y) := minpoly.irreducible hint
  have hne : minpoly F y ≠ 0 := hmonic.ne_zero
  have horb : MulAction.orbit (E ≃ₐ[F] E) y = (minpoly F y).rootSet E := by
    ext z
    rw [Polynomial.mem_rootSet]
    constructor
    · intro hz
      refine ⟨hne, ?_⟩
      have h := (Normal.minpoly_eq_iff_mem_orbit E).mpr hz
      rw [← h]
      exact minpoly.aeval F z
    · rintro ⟨-, hz⟩
      exact (Normal.minpoly_eq_iff_mem_orbit E).mp
        (minpoly.eq_of_irreducible_of_monic hirr hz hmonic).symm
  rw [MulAction.index_stabilizer, horb, Set.ncard_eq_toFinset_card', Set.toFinset_card]
  exact Polynomial.card_rootSet_eq_natDegree (Algebra.IsSeparable.isSeparable F y) (hN.splits y)

def index_stabilizer_eq_natDegree_minpoly.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- `a ∣ b` なら `v_p a ≤ v_p b`。 -/
theorem padicValNat_le_of_dvd {p a b : ℕ} [Fact p.Prime] (ha : a ≠ 0) (hb : b ≠ 0) (h : a ∣ b) :
    padicValNat p a ≤ padicValNat p b := by
  refine (padicValNat_dvd_iff_le hb).mp (dvd_trans ?_ h)
  exact (padicValNat_dvd_iff_le ha).mpr le_rfl

def padicValNat_le_of_dvd.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★★**`k ≥ 1` の段が閉じた**(本ファイルの到達点)。

`E/F` を有限次 Galois、`E` は超距離なノルム体で `Gal(E/F)` が加法的に作用するとする
(★等長性も忠実性も要らない)。`x : E` が `∀ g, ‖g•x − x‖ ≤ Δ` を満たし
`v_p (deg minpoly_F x) = k+1` なら、`y : E` が在って

* `v_p (deg minpoly_F y) ≤ k`(★**深さが 1 段下がる**)、
* `‖y − p • x‖ ≤ Δ`、
* `∀ g, ‖g • y − y‖ ≤ Δ`(★`Δ` は増えない)。

★★体では両辺を `p` で割って `x' := y / p` とすればよく、そのときの損失は
`‖(p : E)‖⁻¹` ちょうどである。★`deg minpoly (y/p) = deg minpoly y`
(`1/p ∈ F` なので安定化群が変わらない)。
★★この定理は `k = 1` も `k ≥ 2` も同じ 1 本で扱う ——
★**持ち場が「`k ≥ 2` は別扱い」と想定していたのは、中間体で降りようとしたためである。** -/
theorem exists_natDegree_minpoly_descent
    {F : Type*} [Field F] {E : Type*} [NormedField E] [IsUltrametricDist E] [Algebra F E]
    [FiniteDimensional F E] [Normal F E] [Algebra.IsSeparable F E]
    {p : ℕ} [Fact p.Prime] {x : E} {Δ : ℝ} (hΔ : 0 ≤ Δ)
    (hx : ∀ g : E ≃ₐ[F] E, ‖g • x - x‖ ≤ Δ) (k : ℕ)
    (hk : padicValNat p (minpoly F x).natDegree = k + 1) :
    ∃ y : E, padicValNat p (minpoly F y).natDegree ≤ k ∧ ‖y - p • x‖ ≤ Δ ∧
      ∀ g : E ≃ₐ[F] E, ‖g • y - y‖ ≤ Δ := by
  have hHidx : (MulAction.stabilizer (E ≃ₐ[F] E) x).index = (minpoly F x).natDegree :=
    index_stabilizer_eq_natDegree_minpoly x
  obtain ⟨Q, y, hQv, hQy, hy1, hy2⟩ :=
    exists_smul_invariant_of_padicValNat_index_succ (p := p)
      (MulAction.stabilizer (E ≃ₐ[F] E) x) (fun h hh => hh) hΔ hx k (by rw [hHidx]; exact hk)
  refine ⟨y, ?_, hy1, hy2⟩
  have hle : Q ≤ MulAction.stabilizer (E ≃ₐ[F] E) y := fun q hq => hQy q hq
  have hdvd : (minpoly F y).natDegree ∣ Q.index := by
    rw [← index_stabilizer_eq_natDegree_minpoly y]
    exact Subgroup.index_dvd_of_le hle
  have hne : (minpoly F y).natDegree ≠ 0 := by
    rw [← index_stabilizer_eq_natDegree_minpoly y]
    exact index_ne_zero_of_finite _
  have := padicValNat_le_of_dvd (p := p) hne (index_ne_zero_of_finite Q) hdvd
  omega

def exists_natDegree_minpoly_descent.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★★**`AxWildDescent` と同じ形まで来た抽象版**(`p` で割った版)。

`exists_natDegree_minpoly_descent` の `y` を `p` で割ると

* `v_p (deg minpoly_F (y/p)) ≤ k`(★`1/p ∈ F` なので安定化群が変わらない)、
* `‖x − y/p‖ ≤ ‖(p : E)‖⁻¹ · Δ`、
* `∀ g, ‖g • (y/p) − y/p‖ ≤ ‖(p : E)‖⁻¹ · Δ`。

★★これは `AxWildDescent K (fun _ => ‖(p : K.closure)‖⁻¹)` の**抽象版そのもの**である
(`Found/PGC/AxTowerDecay.lean:473` と見比べよ)。
★★**残っているのは体の層への移送だけ**で、数学はここで尽きている。
★`‖(p : K.closure)‖⁻¹ = p` なので、これが与える定数は `c k ≡ p`。
★`axDecay p k = p^{(1/(p−1))p^{1−k}}` にするには**分岐**が要る(モジュール docstring の穴 3)。 -/
theorem exists_natDegree_minpoly_descent_div
    {F : Type*} [Field F] {E : Type*} [NormedField E] [IsUltrametricDist E] [Algebra F E]
    [FiniteDimensional F E] [Normal F E] [Algebra.IsSeparable F E]
    {p : ℕ} [Fact p.Prime] {x : E} {Δ : ℝ} (hΔ : 0 ≤ Δ) (hp0 : (p : E) ≠ 0)
    (hx : ∀ g : E ≃ₐ[F] E, ‖g • x - x‖ ≤ Δ) (k : ℕ)
    (hk : padicValNat p (minpoly F x).natDegree = k + 1) :
    ∃ y : E, padicValNat p (minpoly F y).natDegree ≤ k ∧
      ‖x - y‖ ≤ ‖(p : E)‖⁻¹ * Δ ∧ ∀ g : E ≃ₐ[F] E, ‖g • y - y‖ ≤ ‖(p : E)‖⁻¹ * Δ := by
  obtain ⟨z, hz1, hz2, hz3⟩ := exists_natDegree_minpoly_descent hΔ hx k hk
  have hpi : ((p : E))⁻¹ ≠ 0 := inv_ne_zero hp0
  have hsmul : ∀ (g : E ≃ₐ[F] E) (w : E), g • w = g w := fun _ _ => rfl
  have hstab : MulAction.stabilizer (E ≃ₐ[F] E) ((p : E)⁻¹ * z)
      = MulAction.stabilizer (E ≃ₐ[F] E) z := by
    ext g
    simp only [MulAction.mem_stabilizer_iff, hsmul, map_mul, map_inv₀, map_natCast]
    constructor
    · intro h; exact mul_left_cancel₀ hpi h
    · intro h; rw [h]
  refine ⟨(p : E)⁻¹ * z, ?_, ?_, ?_⟩
  · rw [← index_stabilizer_eq_natDegree_minpoly, hstab, index_stabilizer_eq_natDegree_minpoly]
    exact hz1
  · have hxz : x - (p : E)⁻¹ * z = (p : E)⁻¹ * ((p : ℕ) • x - z) := by
      rw [nsmul_eq_mul, mul_sub, ← mul_assoc, inv_mul_cancel₀ hp0, one_mul]
    rw [hxz, norm_mul, norm_inv]
    exact mul_le_mul_of_nonneg_left (by rw [norm_sub_rev]; exact hz2) (by positivity)
  · intro g
    have h : g • ((p : E)⁻¹ * z) - (p : E)⁻¹ * z = (p : E)⁻¹ * (g • z - z) := by
      rw [hsmul, hsmul, map_mul, map_inv₀, map_natCast, mul_sub]
    rw [h, norm_mul, norm_inv]
    exact mul_le_mul_of_nonneg_left (hz3 g) (by positivity)

def exists_natDegree_minpoly_descent_div.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

end ABC3.Found.PGC
