import Mathlib.Analysis.Normed.Order.Lattice
import ABC3.Skeleton.IUTchIII.Cor312

/-!
# [IUTchIII] Corollary 3.12 — 不定性を自明化すると statement は自明になる

`Skeleton/IUTchIII/Cor312.lean` の `cor_3_12` について、
**明らかな退化**を機械にかける。

退化のさせ方は原文の非対称性が指す一点だけ:
Θ 側が q 側と違うのは「**可能な像**の和集合の**正則包**」を測る点なので、

- 可能な像を **`{qImage}` ただ1つ**にする(= 不定性が何も動かさない)
- 正則包を **恒等**にする

とすれば Θ 側と q 側は同じものを測ることになる。

## 結果(下で証明、`sorry` 無し)

`thetaLogVol = qLogVol`(`trivialised_thetaLogVol_eq_qLogVol`)。
したがって Corollary 3.12 の結論は**等号として成り立ち**、
しかも台 `A`・対数体積 `lv`・像 `S`・`|log(q)|` の値が**何であっても**成り立つ
(`trivialised_satisfies_cor_3_12`)。不等式は何も言っていない。

これは Scholze–Stix の指摘

> the critical [IUTT-3, Theorem 3.11] does not become false, but **trivial**

を、**我々のモデルの上で**再現したものである。

## ★これは原典への判定ではない(必読)

自明化したのは我々が置いた `Interface`(`PilotObjectData`)であって、原典ではない。
自明化が起きた原因は次のどれでもありうる:

1. **我々の `Interface` が弱い** — 原文が (Ind1)(Ind2)(Ind3) や正則包に課している
   条件を、我々が逐語から拾い切れていない。`PilotObjectData` は
   `hull` に何の条件も課しておらず、`possibleThetaImages` にも課していない。
2. **必要な数学が未構築** — 条件を書こうにも mono-theta 環境・log-shell・
   Frobenioid・tempered 基本群が無い(mathlib に 0 件)。
3. 原典側の飛躍。

**既定は 1 である**(PLAN §5-2)。3 を名乗るには複数の独立な型設計で同じ壁に当たることと
falsifier が要る。現時点でそれは満たしていないので、`Gap/` には登録しない。

**falsifier(何が起きればこの退化が 1 に落ちるか)**: 原文の逐語から
「`hull` は包含を保つ」「`possibleThetaImages` は空でない」等の条件を
**追加で読み取れた**とき、この退化 witness はそれを満たさなくなりうる。
そのときは本ファイルの witness を作り直す。

**これは原典の主張ではない**(我々のモデルについての事実)ので `.src` を持たない。
-/

namespace ABC3.Check.IUTchIII

open ABC3.Interface.IUTchIII ABC3.Skeleton.IUTchIII

variable {A : Type} (lv : Set A → WithTop ℝ) (S : Set A) (a : ℝ)

/-- **不定性を自明化した `PilotObjectData`**。

`Interface` が課している条件(`qAbs_pos` と `qLogVol_eq`)は残したまま、
Θ 側だけを潰す——可能な像は `{S}` の1つ、正則包は恒等。
台 `A`・対数体積 `lv`・像 `S`・`a = |log(q)|` は**任意**でよい。 -/
noncomputable def trivialised (τ : TopologicalSpace A) (ha : 0 < a)
    (h : lv S = ((-a : ℝ) : WithTop ℝ)) (hS : @IsCompact A τ S)
    (hsub : ∀ U K : Set A, @IsCompact A τ K → U ⊆ K → @IsCompact A τ U)
    (hfin : ∀ U : Set A, @IsCompact A τ U → lv U ≠ ⊤) :
    PilotObjectData where
  Amb := A
  topology := τ
  logVolCompact := fun U h => (lv U).untop (hfin U h)
  logVol := lv
  hull := id
  possibleThetaImages := {S}
  qImage := S
  qAbs := a
  qAbs_pos := ha
  qLogVol_eq := h
  outputLogVolumes := {x : ℝ | (x : WithTop ℝ) ≤ lv S}
  claimsMovedToSkeleton := trivial
  logShellPacket := S
  logVol_eq_of_isCompact := fun U h => (WithTop.coe_untop _ (hfin U h)).symm

/-! ### ★2026-08-15: `Interface` から出した5つの主張を、退化 witness は**全部満たす**

`Skeleton/IUTchIII/Cor312Claims.lean` へ出した主張は、`trivialised` については
すべて成り立つ。すなわち**主張を外に出しても退化 witness は死なない**——
退化を殺すのは主張の置き場所ではなく、`Interface` が対象を同定できているかである。 -/

theorem trivialised_PossibleImagesContained (τ : TopologicalSpace A) (ha : 0 < a)
    (h : lv S = ((-a : ℝ) : WithTop ℝ)) (hS : @IsCompact A τ S)
    (hsub : ∀ U K : Set A, @IsCompact A τ K → U ⊆ K → @IsCompact A τ U)
    (hfin : ∀ U : Set A, @IsCompact A τ U → lv U ≠ ⊤) :
    PossibleImagesContained (trivialised lv S a τ ha h hS hsub hfin) := by
  intro U hU
  have hUS : U = S := hU
  exact hUS.subset

theorem trivialised_LogShellPacketCompact (τ : TopologicalSpace A) (ha : 0 < a)
    (h : lv S = ((-a : ℝ) : WithTop ℝ)) (hS : @IsCompact A τ S)
    (hsub : ∀ U K : Set A, @IsCompact A τ K → U ⊆ K → @IsCompact A τ U)
    (hfin : ∀ U : Set A, @IsCompact A τ U → lv U ≠ ⊤) :
    LogShellPacketCompact (trivialised lv S a τ ha h hS hsub hfin) := hS

theorem trivialised_HullCompactOfRelCompact (τ : TopologicalSpace A) (ha : 0 < a)
    (h : lv S = ((-a : ℝ) : WithTop ℝ)) (hS : @IsCompact A τ S)
    (hsub : ∀ U K : Set A, @IsCompact A τ K → U ⊆ K → @IsCompact A τ U)
    (hfin : ∀ U : Set A, @IsCompact A τ U → lv U ≠ ⊤) :
    HullCompactOfRelCompact (trivialised lv S a τ ha h hS hsub hfin) :=
  fun U K hK hUK => hsub U K hK hUK

theorem trivialised_OutputLogVolumesEq (τ : TopologicalSpace A) (ha : 0 < a)
    (h : lv S = ((-a : ℝ) : WithTop ℝ)) (hS : @IsCompact A τ S)
    (hsub : ∀ U K : Set A, @IsCompact A τ K → U ⊆ K → @IsCompact A τ U)
    (hfin : ∀ U : Set A, @IsCompact A τ U → lv U ≠ ⊤) :
    OutputLogVolumesEq (trivialised lv S a τ ha h hS hsub hfin) := by
  show {x : ℝ | (x : WithTop ℝ) ≤ lv S} = _
  simp [trivialised]

theorem trivialised_QLogVolMem (τ : TopologicalSpace A) (ha : 0 < a)
    (h : lv S = ((-a : ℝ) : WithTop ℝ)) (hS : @IsCompact A τ S)
    (hsub : ∀ U K : Set A, @IsCompact A τ K → U ⊆ K → @IsCompact A τ U)
    (hfin : ∀ U : Set A, @IsCompact A τ U → lv U ≠ ⊤) :
    QLogVolMem (trivialised lv S a τ ha h hS hsub hfin) :=
  le_of_eq h.symm

/-- ★**退化の核心**: 自明化すると Θ 側と q 側は**同じ値**になる。 -/
theorem trivialised_thetaLogVol_eq_qLogVol (τ : TopologicalSpace A) (ha : 0 < a)
    (h : lv S = ((-a : ℝ) : WithTop ℝ)) (hS : @IsCompact A τ S)
    (hsub : ∀ U K : Set A, @IsCompact A τ K → U ⊆ K → @IsCompact A τ U)
    (hfin : ∀ U : Set A, @IsCompact A τ U → lv U ≠ ⊤) :
    thetaLogVol (trivialised lv S a τ ha h hS hsub hfin) = qLogVol (trivialised lv S a τ ha h hS hsub hfin) := by
  simp [thetaLogVol, qLogVol, trivialised]

/-- ★**Corollary 3.12 の結論が、残りのデータが何であっても成り立つ**。

`sorry` 無し。すなわち自明化した `Interface` の下では、
`cor_3_12` は**内容を持たない**。 -/
theorem trivialised_satisfies_cor_3_12 (τ : TopologicalSpace A) (ha : 0 < a)
    (h : lv S = ((-a : ℝ) : WithTop ℝ)) (hS : @IsCompact A τ S)
    (hsub : ∀ U K : Set A, @IsCompact A τ K → U ⊆ K → @IsCompact A τ U)
    (hfin : ∀ U : Set A, @IsCompact A τ U → lv U ≠ ⊤) :
    thetaLogVol (trivialised lv S a τ ha h hS hsub hfin) ≠ ⊤ ∧
      qLogVol (trivialised lv S a τ ha h hS hsub hfin) ≤ thetaLogVol (trivialised lv S a τ ha h hS hsub hfin) := by
  refine ⟨?_, le_of_eq (trivialised_thetaLogVol_eq_qLogVol lv S a τ ha h hS hsub hfin).symm⟩
  rw [trivialised_thetaLogVol_eq_qLogVol lv S a τ ha h hS hsub hfin, qLogVol]
  show lv S ≠ ⊤
  rw [h]
  exact WithTop.coe_ne_top

/-- 「i.e.」の形も自明に成り立つ——`C_Θ ≥ −1` は等号の場合 `C_Θ = −1` で達成される。
すなわち退化した witness は、**abc へ効く形の主張も**空虚に満たす。 -/
theorem trivialised_satisfies_CTheta (τ : TopologicalSpace A) (ha : 0 < a)
    (h : lv S = ((-a : ℝ) : WithTop ℝ)) (hS : @IsCompact A τ S)
    (hsub : ∀ U K : Set A, @IsCompact A τ K → U ⊆ K → @IsCompact A τ U)
    (hfin : ∀ U : Set A, @IsCompact A τ U → lv U ≠ ⊤) (C : ℝ) (hC : thetaLogVol (trivialised lv S a τ ha h hS hsub hfin) ≤ ((C * a : ℝ) : WithTop ℝ)) :
    -1 ≤ C := by
  rw [trivialised_thetaLogVol_eq_qLogVol lv S a τ ha h hS hsub hfin, qLogVol] at hC
  show (-1 : ℝ) ≤ C
  have h4 : -a ≤ C * a := by
    have : ((-a : ℝ) : WithTop ℝ) ≤ ((C * a : ℝ) : WithTop ℝ) := by
      rw [← h]; exact hC
    exact_mod_cast this
  nlinarith [ha]

-- 退化 witness であるにもかかわらず、受理ゲートを全部通ることの確認
#print axioms trivialised_satisfies_cor_3_12

/-! ## ★★ 反証は死んだ — 原文 p.184 と p.175 を足した後の最終測定(2026-08-14)

以前ここには **`cor_3_12` は現在の `Interface` の下では偽である**という定理
(`cor_3_12_refutable_under_current_interface`)と、その反例
(`possibleThetaImages := {∅}`、空集合の対数体積を `⊤` にしたもの)が置いてあった。

原文へ当たり直して**2段階**で足した結果、反例は**構成できなくなった**。

| 追加したもの | 出所 | 効いた先 |
|---|---|---|
| `outputLogVolumes` / `outputLogVolumes_eq` / `qLogVol_mem` | 物理 p.184、Step (xi-e)(xi-f) | **(b) 不等式** |
| `logShellPacket` / `logShellPacket_isCompact` / `possibleThetaImages_subset_logShellPacket` / `hull_isCompact_of_subset_isCompact` / `logVolCompact` + `logVol_eq_of_isCompact` | 物理 p.31 + p.115 + p.127 + p.131 + p.146 | **(a) 有限性** |

結果:

| | 内容 | 結果 |
|---|---|---|
| **(b)** | `qLogVol D ≤ thetaLogVol D` | **証明できる** |
| **(a)** | `thetaLogVol D ≠ ⊤` | **証明できる** |

すなわち `Skeleton.IUTchIII.cor_3_12` は **`sorry` 無しで通る**。

### ★★判定(課題の第3段)— **祝う話ではない**

**(a)(b) の両方が通ったので、我々のモデルの中で Corollary 3.12 は自明になった。**
これは Scholze–Stix の「trivial, not false」を我々のモデル上で再現したものである。

**我々は何も導出していない。** `cor_3_12` の証明は、`Interface` の
`Interface` のフィールドを
そのまま組み合わせただけである。(xi-e)(xi-f) も p.175 のコンパクト性も、
すべて**原文の証明から仮説として輸入**した。

**これは原典への判定ではない。** 原文は実際に証明を書いており、
我々が写していないのはその中身の方である。

### 旧反例はどの条件で死んだか(下で機械的に確定させる)

★`⋃₀ {∅} = ∅` は**コンパクト**なので、コンパクト性の側は**満たしてしまう**。
死因は p.31 の「compact, hence of finite log-volume」に当たる橋
(現在は `logVolCompact` + `logVol_eq_of_isCompact` として型に入っている)の方である。

## ★★★ 1段深い転写の測定(2026-08-14、訂正版)— posit は減ったか

`thetaUnion_isCompact`(= p.175 の「[easily verified] compactness」)は
**原文が検証を書いていない** posit だった。これを1段深く分解し、増減と性質を測った。

### 分解前(2 posit)

| posit | 分類 | 出所 |
|---|---|---|
| `thetaUnion_isCompact` | (β) | p.175。ただし「[easily verified]」で**検証なし** |
| `logVol_hull_ne_top_of_isCompact` | (β) 融合 | p.31 + p.127 を1本に潰したもの |

### 分解後(2 データ + 3 posit)

| posit / データ | 分類 | 出所(すべて 400dpi 目視) |
|---|---|---|
| `logShellPacket : Set Amb` | (β) データ | Rmk 3.9.5, (vii), (Ob1)(p.131)の「tensor packets of log-shells」 |
| `logVolCompact : ∀ U, IsCompact U → ℝ` | (α) データ | **Prop 3.9, (i)(p.115)**: 対数体積は `𝔐(−) → ℝ` |
| `logShellPacket_isCompact` | **(α)** | **p.31 (v) `(a_non)`「I†Fv is compact」** / p.146 |
| `possibleThetaImages_subset_logShellPacket` | **(β)** | **Rmk 3.9.5, (vii), (Ob1)(p.131)** |
| `hull_isCompact_of_subset_isCompact` | **(α)** | **Rmk 3.9.5, (i)(p.127)** の場合分け |
| `logVol_eq_of_isCompact` | (α) 定義的 | `logVol` が `logVolCompact` の拡張であること |

### ★★訂正(2026-08-14)— 一度 (γ) と報告したが誤りだった

`possibleThetaImages_subset_logShellPacket` を「原文に無い(γ)」と報告したが、
**原文にある**。Remark 3.9.5, (vii), (Ob1)(物理 p.131、400dpi 目視):

原文 (IUTchIII p.131):
> various “possible images” that occur as the output of the multiradial al-
> gorithms under consideration are regions — i.e., in essence, elements ∈P
> — contained in tensor packets of log-shells I[scr]_k

しかも (vii) の冒頭は「The operation of forming the hull will play a crucial role
in the context of Corollary 3.12 below」であり、まさに Corollary 3.12 のための箇所である。

**誤りの原因は探索手順**: 我々は Lean 側の名前(`subset` / `contains`)で grep した。
原文は**原文側の語「possible images」**を主語にしていた。
同種の失敗が2回目(1回目は `relatively compact` で探して `compactness` を見落とした)。
規則を `tools/check.mjs` 冒頭 A2 に明文化した。

### 増減のまとめ(訂正後)

- **(α) 0 → 4**(`logVolCompact` / `logShellPacket_isCompact` /
  `hull_isCompact_of_subset_isCompact` / `logVol_eq_of_isCompact`。
  すべて mathlib 語彙で書け、かつ実質的な典拠つき)
- **(β) 2 → 2**(`logShellPacket` データ + `possibleThetaImages_subset_logShellPacket`)
- **★(γ) 0 → 0**

自分で置いた基準「(α) が増えて (γ) が増えないなら収束方向」に照らすと、
**この1段は収束方向**である。

### 副作用: 退化 witness が狭まった

`hull := id` を使う `trivialised` は、新しい `hull_isCompact_of_subset_isCompact` の
ために「コンパクト集合の部分集合はコンパクト」という仮定(`hsub`)を要求するように
なった。分解は**退化の自由度も削った**。

### 有限性は posit ではなく **導出** になった

`logVol_ne_top_of_isCompact` という posit は**廃止**した。Proposition 3.9, (i) が
対数体積を `𝔐(−) → ℝ`(`𝔐(−)` はコンパクト集合)と定めている以上、それは
仮定ではなく**定義域と値域**である。型に写した(`logVolCompact` + `logVol_eq_of_isCompact`)
ので、`Skeleton.IUTchIII.cor_3_12` の有限性の側は
`logVol_eq_of_isCompact` を書き換えて `WithTop.coe_ne_top` で閉じる——**導出**である。
-/

open Classical in
/-- 旧反例が使っていた対数体積:「空集合は `+∞`、それ以外は `−1`」。 -/
noncomputable def lvUnit (s : Set Unit) : WithTop ℝ :=
  if () ∈ s then ((-1 : ℝ) : WithTop ℝ) else ⊤

/-- ★**旧反例の死因を機械的に確定させる**。

`∅` はコンパクトなので `thetaUnion_isCompact` は満たされてしまう。
破れるのは `logVol_hull_ne_top_of_isCompact`(`hull = id` なので `logVol ∅ = ⊤`)。
すなわち反証を殺したのは **p.31 の「compact, hence of finite log-volume」の橋**である。 -/
theorem oldCounterexample_dies_on_the_bridge :
    IsCompact (∅ : Set Unit) ∧ ¬ (∀ U : Set Unit, IsCompact U → lvUnit (id U) ≠ ⊤) := by
  refine ⟨isCompact_empty, fun h => ?_⟩
  exact h ∅ isCompact_empty (by simp [lvUnit])

#print axioms oldCounterexample_dies_on_the_bridge

/-!
## 読み

`cor_3_12` の結論は、自明化した `PilotObjectData` の下で **等号**になる。
不等式 `≥` が主張になるのは、Θ 側の「可能な像」が本当に**複数**あり、
その和集合の正則包が q 側の像より**真に大きい**ときだけである。
その「複数性」を作っているのが (Ind1)(Ind2)(Ind3) だから、
**不定性の内容を書き下さない限り、Corollary 3.12 は内容を持たない。**

これは `Check/PGC/InertiaDegeneracy.lean`(惰性群を自由なデータに置くと
`I_K := ⊥` でも `⊤` でも Corollary 1.3 が成り立つ)と同型の失敗であり、
`Meta/Calibration.lean` が実演した「型は付くが何も言っていない」の3例目である。

**次に何をすべきか**: `possibleThetaImages` と `hull` に、
原文の逐語から**追加で読み取れる条件**があるかを洗う。
無ければ、それは「原文がこの段で条件を与えていない」ことの実測になり、
`Gap/` の候補として §5-2 のトリアージにかける(既定は①、我々のモデル化の誤り)。
-/


/-! ## ★★2026-08-15: 主張を `Interface` から出したとき `sorry` で置いてよいか

「`Interface` にはデータだけを置き、原典が証明している主張は `Skeleton` に
`sorry` で置く」という設計原則を検証する。

結論: **素朴な形は採れない。** `PilotObjectData` はデータの袋であって、
それが原典の対象であることを同定する条件を持っていない。したがって
出した主張は**任意の `D` については偽**であり、`sorry` で置けば
「いつか埋まる借金」ではなく**永久に埋まらない偽の言明**になる。

以下、5つすべてについて**実際に反証を構成する**。意見ではなく測定である。 -/

/-- 反証用の最小 `PilotObjectData`。**データ部分だけ**を埋めてある
(`Interface` に残した整合条件 `qAbs_pos` / `qLogVol_eq` / `logVol_eq_of_isCompact`
はすべて満たす)。台は `ℝ`、対数体積は定数 `−1`。 -/
noncomputable def bareData (imgs : Set (Set ℝ)) (packet : Set ℝ)
    (hl : Set ℝ → Set ℝ) : PilotObjectData where
  Amb := ℝ
  topology := inferInstance
  logVolCompact := fun _ _ => -1
  logVol := fun _ => ((-1 : ℝ) : WithTop ℝ)
  hull := hl
  possibleThetaImages := imgs
  qImage := Set.univ
  qAbs := 1
  qAbs_pos := one_pos
  qLogVol_eq := rfl
  outputLogVolumes := ∅
  claimsMovedToSkeleton := trivial
  logShellPacket := packet
  logVol_eq_of_isCompact := fun _ _ => rfl

/-- ★(Ob1) は反証できる —— 可能な像がパケットに含まれない `D` が作れる。 -/
theorem ob1_is_refutable :
    ¬ PossibleImagesContained (bareData {Set.univ} ∅ id) := by
  intro h
  exact (h Set.univ rfl) (Set.mem_univ (0 : ℝ))

/-- ★パケットのコンパクト性は反証できる —— `ℝ` 全体をパケットに取ればよい。 -/
theorem packetCompact_is_refutable :
    ¬ LogShellPacketCompact (bareData {Set.univ} Set.univ id) := by
  intro h
  exact (NoncompactSpace.noncompact_univ (X := ℝ)) h

/-- ★正則包のコンパクト性は反証できる。 -/
theorem hullCompact_is_refutable :
    ¬ HullCompactOfRelCompact (bareData {Set.univ} ∅ (fun _ => Set.univ)) := by
  intro h
  have he : @IsCompact (bareData {Set.univ} ∅ (fun _ => Set.univ)).Amb
      (bareData {Set.univ} ∅ (fun _ => Set.univ)).topology ∅ := @isCompact_empty ℝ _
  exact (NoncompactSpace.noncompact_univ (X := ℝ)) (h ∅ ∅ he subset_rfl)

/-- ★(xi-e) は反証できる。 -/
theorem outputEq_is_refutable :
    ¬ OutputLogVolumesEq (bareData {Set.univ} ∅ id) := by
  intro h
  have h' : (∅ : Set ℝ)
      = {x : ℝ | (x : WithTop ℝ) ≤ ((-1 : ℝ) : WithTop ℝ)} := h
  have hmem : (-1 : ℝ) ∈ {x : ℝ | (x : WithTop ℝ) ≤ ((-1 : ℝ) : WithTop ℝ)} := by
    simp only [Set.mem_setOf_eq]
    exact le_refl _
  rw [← h'] at hmem
  exact hmem

/-- ★(xi-f) は反証できる。 -/
theorem qMem_is_refutable :
    ¬ QLogVolMem (bareData {Set.univ} ∅ id) := fun h => h

/-! ### 測定結果

出した5つの主張は**すべて反証可能**である。ゆえに:

- `theorem ... := sorry` として置けば、それは**偽の言明**になる。
  posit は「仮定する」と正直に言っているが、偽の `sorry` は嘘である。
- 採った形は**名前付きの `Prop`** + `cor_3_12` の**明示的な仮説**。
  これなら偽を主張せずに、辺の先を**節点**にできる。

★同時に、`trivialised` はこの5つを**全部満たす**
(`trivialised_PossibleImagesContained` 以下)。すなわち
**主張を外に出しても退化 witness は死なない**。退化を殺すのは主張の置き場所ではなく、
`Interface` が対象を同定できているかである。 -/

#print axioms ob1_is_refutable
#print axioms qMem_is_refutable

end ABC3.Check.IUTchIII
