# [CorrHyp] トラックのゴール(2026-09-04 設定、ABC3b)

対象: S. Mochizuki, *Correspondences on Hyperbolic Curves* [CorrHyp]
(物理 18 ページ、著者 1 名、`papers.json` に短縮タグ登記済み、`0_Source` に PDF/txt あり)。

**★2026-09-04、p.3-13・p.15・p.17 を 200dpi 目視確認済み**(`papers.json`
`notationRisk: "medium"`)。実害は overline(X̄ → X)と script 文字(𝒳・𝒴 → X・Y)の
脱落の 2 点(GenEll のような行列の順序入れ替えは無い、詳細は `papers.json` を見よ)。

再現: `node tools/paper-items.mjs CorrHyp`(ABC3c が LocProP のために新設した汎用ツール。
行頭 `Kind N.M.`(ピリオド+行末 or 開き括弧)の宣言規則で番号付き項目を §ごとに数える)。
24 件・重複なしを確認済み(手作業での逐次読解と突き合わせて一致)。

---

## 0. ゴール(現在地)

> **CorrHyp §1 5/5, §2 6/6, §3 3/3, §4 2/2, §5 7/7, §6 1/1 —— 合計 24/24(Skeleton)**

★★**2026-09-04、S0(Skeleton = statement を型で固定する段)を完了した。**
`Interface/CorrHyp/HyperbolicCurve.lean`(`HyperbolicCurveData`・`StackType` を posit)
と `Skeleton/CorrHyp/Section1.lean`〜`Section6.lean` に 24 項目すべてを置き、
`lake build` 0 エラー・`tools/check.mjs` で G1(出典・逐語照合)を全項目パスさせた。
残る `sorry` は 19 件(Definition 5 件は sorry 無しで構成済み)。G9(非空虚性の対照)は
14 件が未着手のまま——プロジェクト全体でトラッキング対象の既知 debt であり、
`check.mjs` 自身が「新規は落とす(ブロックしない)」と明記している。
次の段(S1: Track B で `Found/CorrHyp/` を積む)はまだ手つかず。

## 1. §ごとの内訳(節タイトルつき)

| § | 節タイトル | 件数 | 内訳(宣言順) | 物理ページ |
|---|---|---:|---|---|
| 1 | Basic Definitions | 5 | Definition 1.1 / 1.2 / 1.3 / 1.4 / 1.5 | 3–4 |
| 2 | Review of Results of Margulis and Takeuchi | 6 | Definition 2.1 / 2.2 / 2.3, Proposition 2.4, Theorem 2.5 / 2.6 | 5–7 |
| 3 | The Non-arithmetic Case | 3 | Definition 3.1, Proposition 3.2, Theorem 3.3 | 8 |
| 4 | The Main Theorem | 2 | Lemma 4.1, Theorem 4.2 | 9–10 |
| 5 | Isogenies of General Curves | 7 | Lemma 5.1, Definition 5.2, Theorem 5.3, Lemma 5.4 / 5.5 / 5.6, Theorem 5.7 | 11–15 |
| 6 | Interpretation of a Theorem of Royden | 1 | Theorem 6.1 | 17 |

★`Corollary` は本論文に 0 件。`Remark` は本論文では**すべて無番号**(`"Remark."` のみ、
`Remark N.M` 形式が無い)ので LocProP と同様に項目数に入らない。

★§0(Introduction)は分母 24 から**除外した**(LocProP の §0 とは扱いが異なる——
LocProP の §0 は新規の Definition/Lemma を持つが、本論文の §0 は `Theorem A`/`B`/`C` の
3 件しかなく、そのいずれも本文中で明示的に後続定理の再掲と言明されている。二重計上を
避けるための除外であり、GenEll の §1 起点と同じ扱い):

| §0 の呼称 | 本文中の再掲先の明言(逐語) | ページ |
|---|---|---|
| Theorem A | 「cf. Lemma 4.1 and Theorem 4.2 in the text」 | p.1 |
| Theorem B | 「Theorem B follows from Theorem 5.3 in the text」 | p.2 |
| Theorem C | 「(given as Theorem 6.1 in the text)」 | p.2 |

## 2. 主定理との対応(★本文の言明そのもの、LocProP と違い推定ではない)

| 呼称 | 内容(導入部の記述) | 対応する節 |
|---|---|---|
| Theorem A | 双曲曲線に isogenous な曲線の有限性(第一主定理) | `Lemma 4.1` + `Theorem 4.2`(§4) |
| Theorem B | 一般の曲線は自身の hyperbolic core に一致——非自明な相関を持たない | `Theorem 5.3`(§5) |
| Theorem C | `M_{g,r}` は非自明な自己同型・相関を持たない(Royden の定理の帰結) | `Theorem 6.1`(§6) |

§2(Margulis・Takeuchi の結果の review、`Theorem 2.5`/`2.6`)は他論文 `[Marg]`/`[Take]` の
引用結果であり、本プロジェクトでは証明せず posit する対象になる可能性が高い(要検討——
`Interface` 行き。`genell-track-b` で `Theorem 3.8` 等の外部結果を扱った前例を参照)。

## 3. Track B(Found/CorrHyp/)——2026-09-04 に着手、実測に基づく規模感

★★**Found(sorry 無しの証明)まで24件全部を1セッションで終える見込みは無い**
——GenEll(数か月・現在も§3-§4に3本のsorryが残る)と同型の規模になることを、
実際に mathlib を調べたうえで確認した。「壁」ではなく道として、葉から実装する。

### 使える資産(mathlib、2026-09-04 実測)

| 材料 | 場所 |
|---|---|
| `ℍ`(上半平面)・`SL(2,ℝ)` の Möbius 作用 | `Analysis.Complex.UpperHalfPlane.*` |
| `Subgroup 𝒢 ≤ SL(2,ℝ)`(離散)→ `ℍ` に固有不連続作用 | `UpperHalfPlane.instProperlyDiscontinuousSL2RSubgroup` |
| Schreier の補題(有限指数部分群は有限生成) | `Subgroup.fg_of_index_ne_zero` |
| 被覆空間・モノドロミー | `Topology.Covering.*` |
| `PSL(n,K)` の構成 | `Matrix.ProjectiveSpecialLinearGroup` |

### 無い(§2・§6 の律速)

Margulis/Shimura-arithmeticity(代数群・非可換 Galois コホモロジー・Brauer 群)、
四元数環の分類理論、モジュライスタック `M_{g,r}`、Teichmüller 空間、
Royden の定理——いずれも GenEll のモジュラー多項式ギャップと同型の
「未構築の数学」で、§2(`Definition 2.2`/`2.3`・`Proposition 2.4`・`Theorem 2.5`-`2.6`)
と §6(`Theorem 6.1`)はこれ単体で人年規模になりうる。

### 実装した分(`Found/CorrHyp/FuchsianGroup.lean`、すべて sorry 無し・標準3公理のみ)

★2026-09-04、`Space := FuchsianGroup`(`SL(2,ℝ)` の離散部分群、原文の
`PSL₂(ℝ)⁰` からの逸脱を記録)というモデルで §1・§3 の骨格を実装した:

| 節点 | 実装した宣言 | 状態 |
|---|---|---|
| §1 `Definition 1.1`-`1.4` の土台 | `FEt := IsFiniteIndexIn`・`comp`・`pullback`(`inter`)・`pbFst`・`pbSnd`・`idFEt` | ✅ 圏構造が揃った |
| §1 `Definition 1.5` 直後(isogeny は同値関係) | `isIsogenous_refl`/`_symm`/`_trans` | ✅ **完成** |
| §2 `Definition 2.1` 直前(`Γ ⊆ Comm(Γ)`) | `self_le_commensurator`(mathlib の `Subgroup.Commensurable.commensurator` を発見・流用) | ✅ **完成** |
| §3 `Proposition 3.2`(`Γ_Z ⊆ Comm(Γ_X)`) | `prop_3_2`(`Ncore` = `Z.Γ` の中での `C.Γ` の normal core、を経由) | ✅ **完成**——numbered item の本体として初 |
| §3 `Theorem 3.3` | `fg_of_finiteIndexIn`(Schreier の補題、`Theorem 3.3` が要る片方向) | 部分——下記参照 |

### `Theorem 3.3` に残っていた 2 つの壁

1. **群論側 —— ★2026-09-04 解消**: 原文は「`Comm(Γ_X)` は有限生成部分群 `Γ_X`
   を有限指数に持つので、それ自身も有限生成」と1行で済ませているが、これは
   `fg_of_finiteIndexIn`(有限指数**部分群**が有限生成、Schreier)の**逆向き**
   ——「有限生成な有限指数部分群を持つ群は有限生成」——であり、`exact?` で
   mathlib に見当たらなかった。`fg_of_fg_finiteIndex`(Reidemeister–Schreier の
   最も単純な形、`H` の生成元 + `G/H` の代表系が `G` を生成する)として
   自分で証明した(sorry 無し)。
2. **幾何側(まだ未着手、規模を実測した)**: 「次数が `g',r'` と `e_Y` だけで
   抑えられる」という原文の主張は、`FuchsianGroup` を実際の種数・穴の数
   (§1 冒頭の `type`)に結び付ける Riemann–Hurwitz 型の計算を要求する。
   ★★**mathlib に Gauss–Bonnet・双曲多様体の面積・Riemann 面の種数の理論は
   0 件**(2026-09-04、`grep` で `GaussBonnet`/`RiemannSurface`/幾何的な
   `genus` を実測——`EulerCharacteristic` はホモロジー代数の鎖複体のものだけ)。
   これは GenEll のモジュラー多項式ギャップと同型の「未構築の数学」であり、
   `§3 単体でもTheorem 3.3・Theorem 4.2・Theorem 5.3・Lemma 5.4`-`5.7` 全体の
   共通の律速——**Riemann 面の双曲幾何(面積・種数・Gauss–Bonnet)を
   mathlib レベルで新たに構築する必要がある人年規模の欠落**として記録する。
   `Skeleton/CorrHyp/Section5.lean` の `StackType`/`e_Y` はこの欠落を
   posit で受けている(Skeleton 段階で正しい選択だった)。

### 次の一手(依存の少なさ順)

★★純粋な群論(`FuchsianGroup` の上だけで閉じる節点)は §1・§3 の
`Proposition 3.2` まで**掘り尽くした**——残る節点はどれも「未構築の数学」
(幾何・スキーム論・代数群論)を新たに要求する側に移った:

1. §5 の `StackType`/`e_Y`(Euler 標数)を `FuchsianGroup` の商 `ℍ/Γ` の
   実際の双曲幾何(面積・種数)と結び付ける——mathlib に Gauss–Bonnet が
   無いので、これ自体が新規のライブラリ構築になる。できて初めて
   `Theorem 3.3`/`Theorem 4.2`/`Theorem 5.3`/`Lemma 5.4`-`5.7` に届く。
2. §4(`Lemma 4.1`・`Theorem 4.2`)は係数拡大 `Ext`(k → K)が要り、
   スキーム論(有限型・降下理論)を要求する——①とは別種の未構築の数学。
3. §2 の `MargulisArithmetic`/`ShimuraArithmetic` 本体と §6 の Royden の定理は
   代数群論・Teichmüller 空間という、①②とも異なる第三の未構築の数学を要する。
4. G9(非空虚性の対照)は Track B が本物のデータを供給して初めて
   意味のある形で作れる(`.waiting` に記録済み)。

★どの道も本質的に異なる分野の「未構築の数学」に当たる——GenEll が
モジュラー多項式ギャップ1つで止まったのと違い、本論文は最初から
3つの独立した大きな塊(双曲幾何・スキーム論・代数群論)を要求する
構造になっている(2026-09-04 実測、上記①②③)。

★★**3つの塊のうち②(スキーム論)だけ、mathlib に相応の蓄積がある**
(2026-09-04 追加実測)。`AlgebraicGeometry.Scheme.Etale`・
`LocallyOfFiniteType` に加えて `AlgebraicGeometry/SpreadingOut.lean` に
「有限型なら有限生成部分環まで下りて一意に持ち上がる」という
spreading-out の補題群が存在する。①(双曲幾何、Gauss–Bonnet 0件)・
③(代数群論、非可換 Galois コホモロジー)より見込みがある。

★★★★★★★**2026-09-04、セッション末に①を再実測・確定**
(`genus`・`GaussBonnet`・`RiemannHurwitz` で `grep`——3語とも**文字通り
0件**、無関係な `HurwitzZeta`(数論的ゼータ関数、別分野)のみヒット)。
これで「双曲幾何(種数・Euler標数・Gauss–Bonnet)が mathlib に無い」は
複数回の実測すべてで確認された確定事実——`Theorem 3.3`・§5 の大半・
間接的に `Theorem 4.2`/`Theorem 6.1` の律速。次にここへ戻るセッションは
再確認不要、この人年規模の欠落を埋める作業(mathlib への新規貢献に近い
規模)から直接始めてよい。

★★★**さらに一段実測(2026-09-04)——`Lemma 4.1` に直撃する道具を発見**:
`AlgebraicGeometry/AffineTransitionLimit.lean` は「スキームが affine な
遷移射を持つ余濾過的な図式(`D : I ⥤ Scheme`)の極限である」ときの
spreading-out 定理群を持つ。これは `Spec K`(`K` を有限生成 `k`-部分環
`R_i ⊆ K` の余極限とみなしたときの極限スキーム)の設定そのものであり、
`Lemma 4.1` が要る2つの向きにそれぞれ対応する補題が**両方存在する**:

| `Lemma 4.1` が要ること | 対応する mathlib 補題 |
|---|---|
| `K` 上で定義された対象は、ある有限段 `R_i` まで下りて存在する(spreading out) | `Scheme.exists_π_app_comp_eq_of_locallyOfFinitePresentation` |
| 極限上で一致する2つの射は、ある有限段まで下りればすでに一致する(rigidity/uniqueness) | `Scheme.exists_hom_comp_eq_comp_of_locallyOfFiniteType` |

残る作業は「`X_K`(`K`上の双曲曲線)」「`Z_K` との correspondence」を
この抽象的な極限の枠組みに具体的に載せる段——`I` を `K` の有限生成
`k`-部分環の余濾過圏に取り、`D` をそれぞれの `Spec R_i` 上の底変換の
図式にする、という設定作業がまだ残っている(未着手、規模は未測定だが
道具はすでに mathlib にある)。

★★★★**2026-09-04続報(第23-24件)、上の「設定作業」の骨組みを完成**。
`Found/CorrHyp/FieldLimit.lean`:
- `I := (FgSubalgebra k K)ᵒᵖ`(`K` の有限生成 `k`-部分環全体、包含で有向)
  を実際に構成し、`D := toSchemeDiagram k K := (toRingCat k K).op ⋙ Scheme.Spec`
  という図式を作った。
- **`Spec K` が `D` の極限であること**(`isLimit_specKCone`)を sorry 無しで
  証明した——環側で `K` がその有限生成部分環の余極限であること
  (`isColimitToRingCatCocone`、`Subalgebra.iSupLift` の `RingHom` 版が
  mathlib に無かったので手で構成)を示し、`IsColimit.op` + `Scheme.Spec` が
  極限を保つこと(`Γ ⊣ Spec` の右随伴、`Adjunction.rightAdjoint_preservesLimits`)
  で `Scheme` 側へ運んだ。
- 上の表の2定理が要求する側条件(`IsAffineHom(D.map f)`・
  `CompactSpace(D.obj i)`・`QuasiSeparatedSpace(D.obj i)`)を、
  `D.obj`/`D.map` が常に `Spec(…)` の形であることから**無条件の instance**
  として確立した(mathlib の一般事実だけで済んだ)。

**結果**: `c := specKCone k K`・`hc := isLimit_specKCone k K` として、
`Scheme.exists_hom_comp_eq_comp_of_locallyOfFiniteType`/
`Scheme.exists_π_app_comp_eq_of_locallyOfFinitePresentation` を**側条件抜きで
直接呼べる状態**になった——`Lemma 4.1` のスキーム論的な骨組みは完成。

**まだ残る(誠実な現状)**: `lemma_4_1`(Skeleton)自体は `HyperbolicCurveData`
の抽象フィールド(`Space`・`Ext`・`IsGenericallyScheme`・`ModuliStack`)で
書かれており、これらを「`Spec K` 上の実際のスキーム」に結び付ける
**橋渡しフィールド**が `Interface` にまだ無い。橋渡しを足しても、
`.needs` に記録済みの `LocallyOfFinitePresentation f`(`X_K → Spec k`)や
接空間の計算(標数0での単射性)は依然として要る——モジュライスタック
`M_{g,r}` の実現(§5、mathlib に不在)と分かちがたく結びついている。
現時点で numbered item として Found 完成しているのは
`Proposition 3.2`(§3)のみで、24/24 にはまだ遠い。

★★★★★**2026-09-04続報、`Comm : Fuchsian → Fuchsian` の型付けが
Margulis の二分法そのものを埋め込んでいることを原文照合で確認**(重要な発見)。
`corrhyp_prop_3_2`(`Found/CorrHyp/FuchsianGroup.lean`)に `.src` を付けて
Skeleton の `prop_3_2` と直結させようとして、次に気づいた:

- 原文 p.5(`Comm(Γ)` 導入直後): 「X is not arithmetic. Thus, we have
  Γ ⊆ Comm(Γ) ⊆ PSL2(R)⁰, and **Γ is of finite index in Comm(Γ)**.」
  ——非 arithmetic ⟹ Comm(Γ) 自身が離散(Fuchsian)、という**Margulis の
  二分法の主張そのもの**が §2 の地の文で明言されている。
- `Interface/CorrHyp/HyperbolicCurve.lean` の `Comm : Fuchsian → Fuchsian`
  という型(`Fuchsian` = 離散部分群)は、この事実を**型シグネチャの中に
  既に組み込んでいる**——`Definition 2.1`(`InfinitelyManyCorr`)・
  `Proposition 3.2` の**両方**が `D.Comm Γ : D.Fuchsian`(離散である前提)を
  直接使っているため、これは Skeleton 設計時の見落としではなく、
  原文の記述を素直に型にした結果。
- ★しかし `Comm` が**全域関数**である以上、arithmetic の場合(Comm(Γ) は
  PSL2(R)⁰ で稠密——離散ではない、原文 §2 の対偶側)にも `Fuchsian` 型の
  値を返さねばならず、これは**論理的に無理**——`Comm` は「非arithmeticの
  場合にのみ意味を持つ部分関数」として書き直す必要がある(あるいは
  `MargulisArithmetic`/`ShimuraArithmetic` を経由した場合分けを型に
  持ち込む)。★★これは**Interface の実装可能性そのものに関わる欠陥**
  であり、単なる「橋渡しフィールドが足りない」より深い——`Comm` を
  正しく全域化するには **Margulis の二分法自体を証明する(あるいは
  posit した axiom として明示的に受け取る)** ことが避けられない。
- **`Proposition 3.2` 自身の結論**(`Γ_Z ⊆ Comm(Γ_X)`)は実は
  `Comm(Γ_X)` が離散であることを要求しない(単なる部分群包含の主張)
  ——`corrhyp_prop_3_2` が離散性抜きで証明できたのはこのため。
  ★つまり「`Proposition 3.2` を Skeleton の型のまま Found にする」ことは、
  `Comm` の全域化問題という**§2 由来の別の欠落**を先に解決しない限り、
  型すら整合しない。

**結論(2026-09-04、ユーザーとの対話を経て実施)**: 原文 p.5 の
`Comm(Γ) := {γ ∈ PSL₂(ℝ)⁰ | (γ·Γ·γ⁻¹) ∼ Γ}` を読み直すと、`Comm` の
コドメインを `Fuchsian`(離散)にしたこと自体が誤りで、原文は最初から
「`PSL₂(ℝ)⁰` の部分群」(離散性不問)として定義している。★★**Interface を
実際に修正した**(最小限、加算のみ): `Fuchsian`(離散性を要求しない)は
そのまま残し、`IsDiscrete : Fuchsian → Prop` + `Gamma_isDiscrete : ∀ X,
IsDiscrete (Gamma X)` を新設——`Comm`/`Sub`/`FiniteIndexIn` の型も
`Section2.lean`/`Section3.lean` の記述も一切変えていない(純追加、影響は
ゼロ)。

★★★**さらに一歩進めて、`HyperbolicCurveData` の具体的な項
`corrHypInstance` を `Found/CorrHyp/Instance.lean` に構成した**——
`Space := FuchsianGroup`、`Fuchsian := Subgroup(SL(2,ℝ))`(離散性不問)、
`Gamma`/`Comm`/`Sub`/`FiniteIndexIn`/`IsDiscrete`/`Gamma_isDiscrete` は
`FuchsianGroup.lean` の本物の群論で埋め、`Proposition 3.2` が読まない
残りの field(`MargulisArithmetic` 等)は型合わせの placeholder と明記した。
この `corrHypInstance` において
**`Skeleton.CorrHyp.prop_3_2` の `sorry` を文字通り埋める**
(`prop_3_2_at_instance`、`funext + rfl` で Skeleton の文と関数として
完全一致することを確認済み)ことに成功した——「関連する具体モデルの結果」
ではなく、**Skeleton の主張そのものの実装**という、これまでで最も厳密な
Found 完成(第26件)。

★★★★**2026-09-04続報(第27件)、上の候補を実際に完成させた**。
`Found/CorrHyp/ModularExample.lean`: モジュラー群 `SL(2,ℤ)`(`SL(2,ℝ)` への
像、`discreteSpecialLinearGroupIntRangeSL` で離散)と主合同部分群 `Γ(2)`
(`CongruenceSubgroup.Gamma 2`、mathlib に `FiniteIndex` 済み)を
`corrHypInstance` の `FuchsianGroup` として具体化し、`T=[[1,1],[0,1]]`
witness で両者が**相異なる**ことを確認した上で
`isIsogenous_witness : IsIsogenous corrHypInstance FG_SL2Z FG_Gamma2`
(`sorry` 無し、標準3公理のみ)を得た——`Definition 1.5` が
`corrHypInstance` において**教科書的に意味のある**非空虚性を持つことの
具体的な証拠(有限指数の伝播は `Subgroup.relIndex_map_map` で
`ker(単射) = ⊥` を消し込むだけで閉じた)。

★§2 `Definition 2.1`(`InfinitelyManyCorr`)を `Γ_SL2Z`(モジュラー群)で
非空虚に示せないか検討した(2026-09-04): `SL(2,ℤ)` は古典的に arithmetic
で、その commensurator は `PGL(2,ℚ)`(`Γ_SL2Z` の中で無限指数)——これが
示せれば `Definition 2.1` も非空虚性witnessを得られる。だが必要な
「`diag(2,1/2) ∈ Comm(Γ_SL2Z)`」(Hecke型の共役、`gΓg⁻¹ ∩ Γ` が両側
有限指数)は mathlib に対応する補題が無く(`Hecke`/対角共役で `grep` 0件)、
`Γ(2)` の例のように既存の `FiniteIndex` インスタンスを流用できない——
一から行列計算で `gΓg⁻¹ ∩ Γ` を合同部分群型の集合と同一視する必要がある。

★★★★★**2026-09-04続報(第30件)、完成した**。`g:=diag(2,1/2)`
に対し `g⁻¹Mg = [[a,b/4],[4c,d]]`・`gMg⁻¹ = [[a,4b],[c/4,d]]` を成分計算で
確認し(`conj_inv_formula`/`conj_formula`)、`Γ(4)`(mathlib に
`FiniteIndex` 済み)が両方向の有限指数下界になることを示して
`g2 ∈ Comm(Γ_SL2Z)` を確立、さらに `g2^k`(`k≥1`)の `(1,1)` 成分
`2⁻ᵏ ∈ (0,1)` が整数になりえないことから `g2^k ∉ Γ_SL2Z` を示し、
`{g2^n}` が可算無限個の相異なる剰余類を与えることで
**`Γ_SL2Z` が `Comm(Γ_SL2Z)` の中で無限指数**であることを確立した
(`Found/CorrHyp/HeckeExample.lean`、`infinitelyManyCorr_witness`、
sorry無し・標準3公理のみ)——`Definition 2.1` の非空虚性、§2 の初項目。

★配管の詰め(前回詰まっていた `Matrix.SpecialLinearGroup` の
instances-transparency 問題)は「`set A := (M : Matrix ...)` で台の
`Matrix` を先に1回だけ取り出し、以後 `A` だけで計算する」という
具体的な回避策で解決した——`tools/lean-idioms.md` に追記済み。

★★★★★★**2026-09-04続報(第31-32件)、`Definition 2.3`
(Shimura-arithmetic)を `Γ_SL2Z` について完全に実現した**。
`Found/CorrHyp/ShimuraArithmeticData.lean`: `F:=ℚ`(無限素点1個、
「他で非自明」が空虚に真)・`A:=M_2(ℚ)`(`Algebra.IsCentral.matrix`・
`IsSimpleRing.matrix`・`Module.finrank_matrix` で「4次元中心的単純多元環」
=四元数環という特徴づけを直接満たす、mathlib の基底つき
`QuaternionAlgebra` を経由しない逸脱)・`matrixEquivTensor` で
「唯一の無限素点で自明」・`O_A:=M_2(ℤ)`(`Γ_SL2Z = O_A ∩ SL_2(ℝ)` を
成分レベルで証明)という4データすべてを構成し、`Γ_SL2Z` が
Shimura-arithmetic であることを丸ごと証明した(sorry無し、標準3公理)。
`Definition 2.2`(Margulis-arithmetic)は代数群の部分群スキーム分類・
実点のリー群構造という mathlib に丸ごと不在の理論を要し、引き続き
人年規模。

★次の一手(未着手): `corrHypInstance.ShimuraArithmetic`(現状
`fun _ ↦ False` の placeholder)を `ShimuraArithmeticWitness` へ差し替え
れば、`Definition 2.3` の `.src` を正当に主張できる(道具の集計を
6/24→7/24 に動かせる見込み)。`Instance.lean`(先)→`ModularExample.lean`
→`ShimuraArithmeticData.lean`(後)というimport順のため、差し替えには
ファイル構成の見直しが要る——`prop_3_2_at_instance` 等が
`MargulisArithmetic`/`ShimuraArithmetic` の値を使わないことは確認済み
なので安全なはずだが、未検証。

`Definition 3.1`(`hyperbolicCore`)は今の `core := id` placeholder では
実装にならない(本物の商構成、`Lemma 5.1` の内容)。`Theorem 3.3` は
引き続き Riemann–Hurwitz/Gauss–Bonnet(mathlib 不在)待ち。

★★★**重要な歯止め(2026-09-04 確認)**: `corrHypInstance` の
`MargulisArithmetic`/`ShimuraArithmetic := fun _ ↦ False` という
placeholder を使って `Proposition 2.4`(`MargulisArithmetic Γ ↔
ShimuraArithmetic Γ`)や `Theorem 2.5`(`Arithmetic ↔ InfinitelyManyCorr`)
を「閉じる」ことは**してはいけない**——`Proposition 2.4` は `False ↔ False`
に退化して形式的には閉じるが、これは Margulis/Shimura-arithmeticity の
**実際の内容について何も示していない**(`corrhyp_prop_3_2` のケースと違い、
この2項目は arithmeticity の値そのものが結論に効くので、`prop_3_2` のときの
「未使用の仮定だから placeholder で害はない」という論法が通用しない——
実際 `Theorem 2.5` は my instance では**偽になる**: `SL(2,ℤ)` は現実には
arithmetic で `Comm` の中で無限指数のはずだが、`corrHypInstance` では
`Arithmetic := False` と決め打っているため矛盾する)。★これは
[[report-progress-not-shortfall]] や PGC セッションが記録した「自由な
データによる退化」と同型の罠——`.src` を付けてよいのは、placeholder が
**その項目の結論に一切影響しないことを確認した場合のみ**。

関連: [[leaf-first-with-graph-feedback]] / [[leaves-are-measured-not-guessed]] /
[[measure-mathlib-before-skeleton]] / [[genell-track-b]] / [[corrhyp-track-goal]] /
[[stale-status-read-lean-first]] / [[no-wall-decompose-instead]]

## 4. §4(`Lemma 4.1`・`Theorem 4.2`)の設計案

★§4 は §2/§3/§5/§6 と性質が違う——**mathlib に無い外部定理は不要**
(`isLimit_specKCone`、`Spec K = lim Spec R` が既に sorry 無しで完成しており、
`AffineTransitionLimit.lean` の spreading-out 定理を直接呼べる)。

### 2026-09-04 前半: 関数体案を検討して頓挫(記録として残す)

「`Space := (k 上の関数体) ⊕ (K 上の関数体)`」という和型を検討したが、
`Ext`(係数拡大)がテンソル積 `F ⊗_k K` に対応し、これが**一般には体に
ならない**ことで頓挫した(`Space` 型を全称的に保つのが重くなる)。

### 2026-09-04 後半: Scheme 案に転換、`FEt`/`Ext` を実装完了 ★★★

ユーザーから「工数が多いことで踏みとどまる必要はありません」との明示的な
指示を受け、当初「有限エタール射・ファイバー積の圏論的な整備一式が要る」と
見積もっていた **Scheme 案**に実際に着手したところ、mathlib の基盤が
見積もりよりずっと厚く、**1セッションで完成した**:

- `Found/CorrHyp/SchemeFEt.lean`: `Space := Over BaseK`(`BaseK := Spec ℚ`)。
  `FEt`(有限エタール射)は `AlgebraicGeometry.Etale`/`IsFinite`という
  mathlib の `MorphismProperty` 基盤(合成・恒等・base change 安定性が
  `instance` として自動)にそのまま乗り、`idFEt`/`comp`/`pullback`/`pbFst`/
  `pbSnd` が `infer_instance` と `MorphismProperty.pullback_fst`/`pullback_snd`
  だけで完成。
- `Ext`(係数拡大 `(-)_K`、`K := ℝ`)は、手作りの `Limits.pullback.map` ではなく
  **mathlib 自身の `Over.pullback f ⋙ Over.map f`**(`f : Spec ℝ ⟶ Spec ℚ`)
  という関手の合成として構成すると、`CategoryTheory.MorphismProperty.
  overPullbackMap` が base change 安定性からの遺伝を一発で与えてくれた
  ——`extFEt`(`Ext` の射側の作用)がこれで完成。
- `HyperbolicCurveData.Space : Type` を `Type u`(universe 多相)に変更
  (逸脱、`HyperbolicCurve.lean` に記録)——`Over BaseK : Type 1` が
  旧来の `Type`(`= Type 0`)固定と衝突したため。既存の `corrHypInstance`/
  `corrHypInstance2`(`u=0`)には影響なし。
- `Found/CorrHyp/Instance3.lean`: 上記を全部使い `corrHypInstance3 :
  HyperbolicCurveData` を具体的に構成(`Lemma 4.1` が読まないフィールドは
  安全な placeholder)。`Lemma 4.1` の statement が `corrHypInstance3` に
  対して型検査を通ることを確認。

**結論: `Space`/`FEt`/`Ext` の一貫した具体化という設計上のブロッカーは
解消した。** 残るのは `Lemma 4.1` **本体**(spreading-out の実際の証明)。

### `Lemma 4.1` 本体に向けた実測(2026-09-04、続き)

`toSchemeDiagram ℚ ℝ : (FgSubalgebra ℚ ℝ)ᵒᵖ ⥤ Scheme`(`FieldLimit.lean`)は
`AffineTransitionLimit.lean` が要求する前提を**両方とも自動で満たす**ことを
確認した:
```
example : IsCofiltered (FgSubalgebra ℚ ℝ)ᵒᵖ := inferInstance   -- OK
example {i j} (f : i ⟶ j) : IsAffineHom ((toSchemeDiagram ℚ ℝ).map f) := by
  infer_instance   -- OK
```
`AffineTransitionLimit.lean` には「極限の中で2つの射が一致する⟹有限段階で
一致する」型の spreading(`Scheme.exists_hom_hom_comp_eq_comp_of_
locallyOfFiniteType` 等)・「アフィン開被覆は有限段階に降りる」
(`Scheme.exists_isOpenCover_and_isAffine`)・「アフィン性は有限段階に降りる」
(`Scheme.exists_isAffine_of_isLimit`)は**すでにある**。

★★**まだ無いもの**(`lemma_4_1.needs` の1番目、EGA IV 8.8-8.10 相当の核心):
「極限上の**新しいスキーム**(有限型・有限エタールという条件つき)が、
ある有限段階の**新しいスキーム**の base change として存在する」という
**構成的な**降下——上の3つはどれも「すでにある対象・射」の spreading で、
「新しい対象を`Z_K`から`Z`へ**作る**」ものではない。組み立て方の見通し:
`exists_isOpenCover_and_isAffine` でアフィン開被覆を有限段階に降ろし→
各アフィン片を `RingHom.EssFiniteType.exists_eq_comp_ι_app_of_isColimit`
(`Algebra/Category/Ring/FinitePresentation.lean`)型の環準同型の spreading で
有限段階の環に降ろし→貼り合わせデータ(遷移射の一致)を
`exists_hom_hom_comp_eq_comp_of_locallyOfFiniteType` で有限段階に降ろし→
`Etale`/`IsFinite` の条件自体も有限段階で成り立つことを言う(これは
`RingHom.Etale`/`FormallyEtale` の「有限生成部分環での近似」相当、
未確認)。**これは1つの既存補題では閉じず、組み立てが要る**——
`lemma_4_1.needs` の2番目(étale 剛性)・3番目(標数0接空間)は
まだ手つかず。

★★★2026-09-04さらに続報: `RingHom.Etale`(`FormallyUnramified ∧ Flat ∧
FinitePresentation` に分解可能、`Algebra.Etale.iff_flat_and_
formallyUnramified`)が余極限に沿って spreading する**汎用**定理は
mathlib のインデックスに見当たらない——`FinitePresentation` 部分は
`EssFiniteType.exists_eq_comp_ι_app_of_isColimit` を使い回せるが、
`Flat`/`FormallyUnramified` 部分は個別に確認が要る(未着手)。

★★★★**設計上の発見(2026-09-04、`corrHypInstance3` で `lemma_4_1` を
直接 `intro` して発覚)**: `Skeleton/CorrHyp/Section4.lean` の `lemma_4_1` は
結論を `∃ Z, ∃ (h : ZK = D.Ext Z), ...`(**命題的な等号 `=`**)で書いている。
一方 `HyperbolicCurveData.Iso : Space → Space → Prop` は witness(実際の
同型射)を一切持たない**裸の Prop**。`FuchsianGroup`(`Subgroup`、
extensional な等号を持つ)モデルではこの `=` は無害だったが、**`Scheme`
モデルでは一般に強すぎる**——`ZK` が `Ext X` と同型な**別の構成**
(例: 環を `ULift`/`Spec` で作り直したもの)であっても、finite étale な
correspondence 自体は(同型射も finite étale なので)問題なく存在しうる。
このとき結論の「`ZK = D.Ext Z` を満たす `Z` が存在する」は、`ZK` が
`Ext(-)` の形に**命題的に**一致する保証が無い限り**偽になりうる**——
これは EGA IV の降下理論を組んでも解決しない、**独立した形式化上のギャップ**
(数学的な困難ではなく、`Iso` が witness を持たないという Interface 設計の
issue)。

**次への示唆**: この issue の解消(`Iso` に実際の同型射データを持たせる、
または `lemma_4_1` の結論を `D.Iso` ベースに書き換える)は `Lemma 4.1` の
**構成**そのものを不要にはしない(`Z` を実際に作る EGA IV の仕事は
`=` でも `Iso` でも変わらず必要)——ただし構成できた後の「一致」の
証明が `rfl`/構成的等号ではなく同型の構成で済むようになる、という点で
後工程を楽にする。`Interface/CorrHyp/HyperbolicCurve.lean` を触る前に、
他の numbered item(`FinitelyManyUpTo` 経由で `Iso` を使う `Theorem 3.3`・
`Theorem 4.2`・`Theorem 2.6`)への影響も見る必要がある——`Iso` を
Prop のままにするか、witness 付きの構造体に変えるかは要検討。

次にここへ戻るときの出発点: 上の「構成的な降下」を`FEtK`(有限エタール、
`IsFinite ∧ Etale`)の場合に特化した1つの補題として切り出し、
`exists_isOpenCover_and_isAffine` から着手する。

### 2026-09-04さらに続報: `StandardEtalePair` の降下が完成 ★★★

上の「組み立て方の見通し」の最初の環——標準エタール表示(`f`・`g` の2多項式・
`cond`)の降下——を実際に `FieldLimit.lean` に実装し、sorry 無しで完成させた:

- `exists_fg_subalgebra_polynomial`(1変数)・`_pair`(2変数、同じ部分環に
  同時に降ろす)・`_pair_monic`(monic 性も遺伝)・`_family`(有限個の多項式の
  一般化)——多項式は係数が有限個なので、ある有限生成 `k`-部分環へ必ず
  降ろせるという初等的だが必要な部品。
- `exists_fg_subalgebra_standardEtaleCond`——`cond`(`f'*p₁+f*p₂=g^n`)自体は
  `f, g, p₁, p₂` を`_family`で同時に降ろした上で、`Polynomial.map`
  (`algebraMap R K` が単射なので単射)の単射性で等式を引き戻す。
- `exists_fg_subalgebra_standardEtalePair`——上の2つを束ね、
  `Algebra.StandardEtalePair K` 全体がある有限生成 `k`-部分環上の
  `StandardEtalePair` の像であることを示す。

**残る部品**(`Algebra.StandardEtalePresentation` まで完成させるには):
(a) `x : S`(`lift` が bijective になる元)自体の降下——`S` 自体が
「有限生成部分環上のある `S₀` の base change」であることも言う必要があり、
これは `S`(有限エタール多元環)を `R[X]/(f)` のこの `x` による同型で
"再定義" することで回避できる可能性がある(`x` を明示的に構成する必要が
無くなる)、(b) 一般の有限エタール多元環は「至る所で」標準エタールとは
限らず(`Algebra.basicOpen_subset_etaleLocus_iff_etale`)、**局所的に**
しか標準エタール表示を持たない——`Z_K` 全体のアフィン開被覆
(`exists_isOpenCover_and_isAffine`)の**各アフィン片**を、さらに
étale-locus のレベルで細分してから、この降下を1片ずつ適用し、最後に
`RingHom.EssFiniteType.exists_comp_map_eq_of_isColimit`(遷移射の一致の
降下)で貼り合わせる、という二重の被覆の組み立てが要る。

### 2026-09-04さらに続報: (a) の見通しを精緻化——「テンソル積の余極限保存」は不要

`exists_fg_subalgebra_standardEtalePair` で得た `f₀, g₀ : Polynomial R.1`
(`f₀.map = f`・`g₀.map = g`)から、**`S₀ := (R.1[X]/(f₀))` を `g₀` で局所化
したもの**を直接定義すればよい——これが「有限生成部分環上のある `S₀`」の
候補で、`x` を明示的に探す必要は無い(原文の `x`・`hasMap`・`lift_bijective`
は「`S` が既にこの形をしている」ことの確認用データであり、**こちらから
`S₀` を作る側では不要**)。

残る唯一の作業は「`S₀ ⊗_R K ≅ S`」——これは
`(R[X]/(f₀)) ⊗_R K ≅ K[X]/(f₀.map) = K[X]/(f)`(商環の base change、
`Ideal.map`/`Polynomial.quotientEquiv` 相当)と、局所化の base change
(`IsLocalization.baseChangeEquiv` 相当)を合成すればよい、**標準的な
可換環論の事実**——これは一般の「余極限の保存」(`Under.pushout` が
left adjoint であることを使う圏論的な議論、`Adjunction.leftAdjoint_
preservesColimits`)よりずっと軽い道具で済む可能性が高い。

★`Under.pushout`(`CategoryTheory.Under.pushoutIsLeftAdjoint`、
`PreservesColimit` が `infer_instance` で自動的に付くことは確認済み)による
「テンソル積は左随伴だから余極限を保存する」という一般論も**筋は良い**が、
`toRingCat k K` を `Under (CommRingCat.of k)` へ持ち上げる段で
`Under.mk`/`Under.homMk` の defeq 周りの配管が詰まった(`Over.mk` で
見た「`instances` 透明度」問題の `Under` 版、未解決のまま持ち越し)——
`corrhyp-goal.md` に記録だけして、上記のより軽い道(商環・局所化の base
change)を先に試す方が近道と判断した。

### 2026-09-04さらに続報: `StandardEtalePair.Ring` の base change が完成 ★★★★

上記の見通しどおり、`S₀ := R[X]/(f)` を局所化する代わりに mathlib 自身の
`StandardEtalePair.Ring`(`Polynomial (Polynomial R) ⧸ span {C f, X*C g-1}`、
明示的な二変数多項式商)を直接使う道で、base change の同型を
**完全に(sorry 無しで)構成した**:

- `Bivariate_equivMvPolynomial_map`——`Polynomial.Bivariate.equivMvPolynomial`
  (二変数多項式環 `≃ MvPolynomial (Fin 2)`)が係数の写像 `φ` と可換である
  ことを、両辺が環準同型であることから `Polynomial.ringHom_ext'` を2回
  (外側の変数・内側の変数)適用して示した。
- `standardEtalePair_ring_baseChange`——`equivMvPolynomialQuotient` が使う
  イデアル(`{C f, X*C g-1}` の span)の生成元の像が、`includeRight` →
  `MvPolynomial.algebraTensorAlgEquiv`(多変数多項式環の base change、
  mathlib 既存)で運ぶと、対応する `K` 側の生成元に**文字通り一致する**
  ことを `Ideal.map_span` + 上の可換性で計算した。
- `standardEtalePairRingBaseChange`——上の2つと
  `Algebra.TensorProduct.tensorQuotientEquiv`(商環の base change)・
  `Ideal.quotientEquivAlg`(イデアルの対応から商の同型を作る)を合成し、

  **`TensorProduct R K P₀.Ring ≃ₐ[K] (P₀.map (algebraMap R K)).Ring`**

  を得た。「余極限の保存」という一般論は不要で、mathlib の base change
  補題(`MvPolynomial.algebraTensorAlgEquiv`・`tensorQuotientEquiv`・
  `Ideal.quotientEquivAlg`)の組み合わせだけで閉じた。

**これで `Lemma 4.1` の構成的降下のうち、最も技術的に重かった「1つの
標準エタール表示の base change」が完全に閉じた。** 残るのは組み立てのみ:
(a) 一般の有限エタール多元環は「至る所で」標準エタールとは限らないため
(`Algebra.basicOpen_subset_etaleLocus_iff_etale`)、`Z_K` 全体のアフィン
開被覆(`Scheme.exists_isOpenCover_and_isAffine`)の各片をさらに
étale-locus のレベルで細分してから、この降下を1片ずつ適用すること、
(b) 各片の降下結果を `Scheme.exists_hom_hom_comp_eq_comp_of_
locallyOfFiniteType`(遷移射の一致の降下)で貼り合わせて `Z` 全体を
構成すること。数学的な核心(base change の構成)は完成したので、残るのは
**スキーム論の定型的な貼り合わせ作業**という段階に到達した。

### 2026-09-04さらに続報: `Spec K = lim Spec R` を `Over BaseK` へ持ち上げ完成

`Found/CorrHyp/ExtLimit.lean`(新規): `isLimit_specKCone`(裸の `Scheme` での
極限)を **`Over BaseK`**(`k`-スキーム全体の圏)へ持ち上げた
(`isLimit_specKConeOver`)。`toSchemeDiagramOver ⋙ Over.forget BaseK =
toSchemeDiagram ℚ ℝ` が `rfl` で成り立つ(持ち上げが「余計な情報を足して
いない」ことの直接確認)ことを利用し、`isLimit_specKCone` 自身の
`lift`/`fac`/`uniq` から `Over` 版の `IsLimit` を直接構成した(一般の
「`Over.forget` は連結図式の極限を creates する」という instance
(`CostructuredArrow.Over.createsLimitsOfShapeForgetOfIsConnected` 等)を
`infer_instance` で探したが見つからず、直接構成に切り替えた——これも
「一般論より具体構成」という今回のセッションの教訓の再確認)。

**残る唯一の詰まり**: `Over.pullback X.hom`(`X.hom` に沿った base change、
`Over.mapPullbackAdj` の右随伴なので `PreservesLimit` が `infer_instance`
で自動的に付くことは確認済み)を `specKConeOver` に適用すれば
`Ext X` の極限表示が得られるはずだが、`Limits.pullback X.hom toBaseK`
(`Ext X` の定義で使う引数順)と `Limits.pullback toBaseK X.hom`
(`Over.pullback X.hom` が自然に与える引数順)が入れ替わっている——
`CategoryTheory.Limits.pullbackSymmetry`(mathlib既存、`pullback f g ≅
pullback g f`)で埋める見込み。埋まれば、`FieldLimit.lean` の
`standardEtalePairRingBaseChange`(環側の base change)と合わせて、
`Lemma 4.1` の構成的降下に必要な**スキーム側・環側の両方の道具が揃う**。

### 2026-09-04さらに続報: `isLimit_extCone` 完成 ★★★★★——スキーム側の道具が完全に揃った

`pullbackSymmetry` を経由する代わりに、`Ext X` の定義そのものと**同じ
引数順**(`pullback X.hom (D R).hom`)で図式 `extDiagram X`(`Over.pullback`
を経由しない、`Limits.pullback.map` を直接使う手作りの図式)を組み直し、
その頂点 `Limits.pullback X.hom toBaseK`(= `(Ext X).left`)が極限である
ことを**完全に(sorry 無しで)証明した**(`isLimit_extCone`):

- `extConePi`(頂点への `π`)・`extCone`(cone 本体)・`auxCone3`
  (`Cone.postcompose` 経由で `toSchemeDiagram ℚ ℝ` の cone へ変換)・
  `R0hom`(`⊥` を代表元に取るための射)・`extCone_fst_const`(`X.left` への
  射影が代表元に依らず一定)という6個の補題を積み上げた上で、
  `IsLimit.mk` の `lift`/`fac`/`uniq` を直接構成。
- 配管上の最終ブレークスルーは2つ: (1) `Cone` の `π` フィールドには
  **独立に作った `NatTrans`(`extConePi`)をそのまま代入**し、`naturality`
  を再証明しない、(2) `IsLimit.mk` の `uniq` が与える仮定 `hm` を
  `have hm' : ∀ R, m ≫ (extConePi X).app R = s.π.app R := hm` のように
  **型を明示して再束縛**することで `rw` の構文一致ではなく `have` 自体の
  defeq チェックを経由させる——これで `Functor.const` の配管(`tools/
  lean-idioms.md` 第22項)が完全解決した。

**これで `Lemma 4.1` の構成的降下に必要な道具が完全に揃った**:
`standardEtalePairRingBaseChange`(環側の base change、`FieldLimit.lean`)
+ `isLimit_extCone`(スキーム側の極限表示、`ExtLimit.lean`)。残るのは
純粋な組み立て作業のみ: (a) 一般の有限エタール多元環は「至る所で」標準
エタールとは限らないため、`Z_K` のアフィン開被覆を étale-locus のレベルで
細分してから各片に `exists_fg_subalgebra_standardEtalePair_map` を適用する
こと、(b) 各片の降下結果を `Scheme.exists_hom_hom_comp_eq_comp_of_
locallyOfFiniteType`(遷移射の一致の降下)で貼り合わせて `Z` 全体を構成
すること。数学的な核心はすべて完成した——残るのはスキーム論の定型的な
貼り合わせという最終段階。

### 2026-09-04さらに続報: `AffineTransitionLimit.lean` 適用に要る追加の前提を発見

上の (a)/(b)(アフィン開被覆の細分・遷移射の貼り合わせ)は
`AffineTransitionLimit.lean` の `Scheme.exists_isOpenCover_and_isAffine`
等を `isLimit_extCone`(`extDiagram X` の極限)に**そのまま**適用すれば
済むと想定していたが、これらの補題は `[∀ {i j} (f : i ⟶ j), IsAffineHom
(D.map f)]`(図式の遷移射がアフィン)・`[∀ i, CompactSpace (D.obj i)]`・
`[∀ i, QuasiSeparatedSpace (D.obj i)]` を要求する。`toSchemeDiagram ℚ ℝ`
自身については(`Spec R`が常にアフィン・準コンパクトなので)これらは
`FieldLimit.lean` で既に確立済みだが、**`extDiagram X`(`X.hom` に沿って
base change した図式)については未確認**——`(extDiagram X).map h`
(`Limits.pullback.map` で組んだ射)は `infer_instance` では
`IsAffineHom` が自動的に付かない(`SchemeFEt.lean` の `Ext`/`extFEt` で
`pullback.map` ではなく `Over.pullback`+`overPullbackMap` に乗り換えて
解決したのと同じ理由——`pullback.map` は `MorphismProperty.pullback_fst`/
`pullback_snd` の形に直接一致しない)。

数学的には真(アフィン射の base change はアフィン、という標準事実)なので
**証明できるはずだが未着手**——`(extDiagram X).map h` を `pullback.fst`/
`pullback.snd` の合成として書き直すか、`extDiagram` 自体を
`Over.pullback`(`X.hom` 側)ベースで組み直す(`Ext`/`extFEt` で使った
設計)必要がある。さらに `CompactSpace`/`QuasiSeparatedSpace` は `X.left`
自身にこれらの性質を要求する可能性があり(`X` が一般の `Over BaseK` の
元では保証されない)、`corrHypInstance3` の `Space` を「有限型
(qcqs)なスキームに限る」よう絞る設計変更が要るかもしれない——双曲曲線は
常に有限型なので原文の意図には合致するが、`Interface` 側の変更が要る。
次にここへ戻るときの具体的な出発点。

### 2026-09-04さらに続報: `IsAffineHom` は解決、`qcqs` 前提の必要性は確定

`extDiagram_map_isAffineHom`(`ExtLimit.lean`)として、`extDiagram X` の
遷移射がアフィンであることを**sorry 無しで証明した**——
`CategoryTheory.MorphismProperty.pullbackMap`(`P i₁`・`P i₂` から
`P (pullback.map f g f' g' i₁ i₂ (𝟙 _) ...)` を直接与える、`Over.pullback`
+`overPullbackMap` より単純な道)に `IsAffineHom`(`isAffineHom_
isStableUnderBaseChange` が mathlib 既存)を適用するだけで閉じた。

一方 `CompactSpace ((extDiagram X).obj R)` は `infer_instance` で自動的に
付かないことを直接確認——`(D R).hom`(`Spec R → BaseK`)はアフィンだが
`X.hom` 自身には何の制約も無いため、**`X.left` 自身の準コンパクト性
(=有限型性の一部)を要求する**ことが確定した。これは数学的な困難ではなく
**設計上の要求**(双曲曲線は原文の定義から常に有限型なので、`Lemma 4.1`
を正しく述べるには本来 `X` にこの仮定が必要——mathlibが要求しているのは
paper自身が暗黙に仮定している事実そのもの)。次にここへ戻るときは、
`corrHypInstance3` の `Space` を qcqs スキームの部分型に絞るか、
`lemma_4_1` を適用する具体的な文脈で `X` に `CompactSpace X.left`/
`QuasiSeparatedSpace X.left` を仮定として追加する、のどちらかを選ぶ
設計判断から始める。

### 2026-09-04さらに続報: `QcqsSpace`・`corrHypInstance4`・被覆の降下が完成

`Space` を qcqs スキームの部分型に絞る側を選び、`QcqsSpace.lean`
(`FEt`/`Ext` の圏構造を `QcqsSpace` 上に再構築、すべて sorry 無し)・
`Instance4.lean`(`corrHypInstance4`、`Lemma 4.1` の statement が型検査を
通ることを確認)を完成させた。

さらに `Lemma 4.1` の構成的降下で最も中心的な一手——**`Ext X` の有限アフィン
開被覆が、ある有限段階 `X.left ×_{BaseK} Spec R` の有限アフィン開被覆から
来ること**——を `exists_extDiagram_finite_affine_descent`(`ExtLimit.lean`)
として完全に証明した。鍵となった `Scheme.exists_finite_affineOpenCover`
(任意の準コンパクトなスキームが有限アフィン開被覆を持つという、CorrHyp に
依存しない一般的な事実、`Scheme.affineCover` + `IsCompact.elim_finite_
subcover` で構成)は、`exact?`/`apply?` では見つからず自分で組み立てた。

**残る組み立て**: (a) `exists_extDiagram_finite_affine_descent` で得た
各アフイン片(`X.left ×_{BaseK} Spec R` のアフィン開)上の有限エタール
多元環を étale-locus のレベルでさらに細分し、`exists_fg_subalgebra_
standardEtalePair_map`/`standardEtalePairRingBaseChange`(`FieldLimit.lean`)
で有限段階へ降下すること、(b) 各片の降下結果を `Scheme.exists_hom_hom_
comp_eq_comp_of_locallyOfFiniteType`(遷移射の一致の降下)で貼り合わせて
`Z` 全体を構成すること。`Lemma 4.1` の証明に使う道具はほぼすべて揃った
——残るのは「一般の有限エタール射を局所的に標準エタールへ帰着する」という
最後の理論的なピースと、その後の純粋な貼り合わせ作業。

### 2026-09-04さらに続報: étale-locus 細分(最後の理論的ピース)も完成 ★★★★★

上で「最後の理論的なピース」と記した (a) の étale-locus 細分——
「一般の有限エタール多元環は至る所で標準エタールとは限らない」という
ギャップ——を `exists_finite_standardEtaleCover`(`FieldLimit.lean`)として
完全に解決した:

**étale な多元環は、有限個の基本開集合の上で標準エタールになる。**

`Algebra.IsEtaleAt.exists_isStandardEtale`(mathlib既存、各点で局所的に
標準エタールになる)を `PrimeSpectrum S` の各点に適用し、その準コンパクト性
(`PrimeSpectrum.compactSpace`、任意の環に対して常に成り立つ)で有限部分
被覆に絞り込む——`Scheme.exists_finite_affineOpenCover` と全く同じ
「点ごとの局所性→準コンパクト性で有限化」というパターンが、スキームの
開被覆でも環の基本開被覆でも同様に機能することを確認した。CorrHyp に
一切依存しない一般的な事実。

**これで `Lemma 4.1` の局所的な構成的降下(1つのアフィン片について)に
使う道具がすべて揃った**: `exists_extDiagram_finite_affine_descent`
(被覆の降下)→ `exists_finite_standardEtaleCover`(étale-locus細分)→
`exists_fg_subalgebra_standardEtalePair_map`+`standardEtalePairRingBase
Change`(各標準エタール片の環レベル降下)。残るのは、これらを実際に
組み合わせて1つのアフィン片について局所的な `Z` の構成を完成させること、
そして複数の片を `Scheme.exists_hom_hom_comp_eq_comp_of_
locallyOfFiniteType`(遷移射の一致の降下)で貼り合わせることの2つ——
どちらも「新しい道具を探す」段階ではなく「既存の道具を組み合わせる」
段階に入った。

### 2026-09-04さらに続報: 貼り合わせ機構(`Scheme.GlueData`)の在庫を確認

`Z` を局所片から貼り合わせる最終段で使う道具の在庫を確認した——mathlib は
`AlgebraicGeometry/Gluing.lean` に `Scheme.GlueData`(貼り合わせデータの
構造体)・`Scheme.Cover.glueMorphisms`(局所的に定義された射を1つの
射へ貼り合わせる)・`Scheme.Cover.gluedCover`/`gluedScheme` という
一通りの貼り合わせ機構をすでに持っている——ゼロから構築する必要はない。

**§4(`Lemma 4.1`)の現状のまとめ**: 数学的に必要な道具(環側の base
change、スキーム側の極限表示・被覆の降下・étale-locus 細分、貼り合わせ
機構)は mathlib 既存分・このセッションで新規に証明した分を合わせて
**すべて確認できた**。残るのはこれらを実際に1本の証明として組み立てる
エンジニアリング作業のみ——「必要な数学が見つからない」という壁は
存在しない。組み立て自体の分量はなお大きく、複数の作業単位(1つのアフィン
片の局所構成→複数片の整合性→貼り合わせ→`Z` の全体像の確認)に分解して
段階的に進める必要がある。

### 2026-09-04さらに続報: 「スキームレベル→環レベル」の最後の橋渡しが完成

`ExtLimit.lean` に `Etale.algebraEtale_appLE` を追加(★sorry 無し・標準3
公理のみ)。`[Etale α] [IsFinite α]` なスキーム射 `α : C ⟶ Y` をアフィン開
`U ⊆ Y` へ制限すると、誘導される環準同型 `α.appLE U (α ⁻¹ᵁ U) le_rfl` の
`toAlgebra` により `Algebra.Etale Γ(Y,U) Γ(C, α⁻¹U)` が成り立つ、という
補題——`α ⁻¹ᵁ U` のアフィン性は `IsAffineHom.isAffine_preimage` から
(`[IsFinite α]` → `IsAffineHom`)。

つまずいた点: `Etale.etale_appLE` の実際のシグネチャは
`∀ {X Y} (f : X ⟶ Y) [Etale f] {U} (hU : IsAffineOpen U) {V} (hV : IsAffineOpen V)
(e : V ≤ f⁻¹ᵁU), (f.appLE U V e).hom.Etale` で、射 `f` 自身が(instance
引数より前に)**明示引数**になっている。これを渡さず `Etale.etale_appLE
hU hV le_rfl` と書くと `hU`/`hV`/`le_rfl` がそれぞれ 1 つずつズレた
引数スロットに入ってしまい意味不明な型エラーになる——`Etale.etale_appLE
α hU hV le_rfl` と `α` を明示すれば解決(`tools/lean-idioms.md` に追記
予定)。また `(f.appLE U V e).hom.Etale`(`RingHom.Etale`)は
`Algebra.Etale`(`f.hom.toAlgebra` の下)と単なる `def` の展開で defeq
なので、`letI` で同じ `toAlgebra` インスタンスを立てておけば `exact` が
そのまま通る——追加の変換補題は不要だった。

これで「相関 `c.α`(スキームレベルの有限エタール射)→ アフィン片への
制限 → 環レベルの `Algebra.Etale`」という、`exists_finite_standardEtaleCover`
を実際に呼び出すために必要だった最後の理論的ギャップが埋まった。
`Lemma 4.1` の証明を組み立てるのに必要な道具は、これで**文字通りすべて
Lean の宣言として存在する**状態になった。

### 2026-09-04さらに続報: 実際の組み立てに着手して見つかった2つの発見

`Lemma 4.1` を `corrHypInstance4` に対して実際に組み立てようとして、
statement を精査した結果、2つの重要な発見があった——**どちらも
「見つからない数学」ではなく「statement 自体の精査から出た発見」**
という点で、これまでの発見(EGA IV道具の在庫確認)とは性質が異なる。

**発見1(要修正、未修正): `Definition 1.1`(`Corr`)から原文の
「C is nonempty」が脱落している。** 原文 (p.3) は「a correspondence
...where we assume that C is nonempty」だが、`Skeleton/CorrHyp/Section1.lean`
の `structure Corr` にはこの仮定が無い(このセッションで見つけたが、
`Section1.lean` 自体は以前のセッションで書かれたもの)。これは
**単なる見た目の欠落ではなく実害がある**——mathlib で確認したところ、
空スキーム `∅ = ⊥_ Scheme` から任意の `Y` への唯一の射
`initial.to Y` は `IsFinite.instOfIsClosedImmersion`(閉埋め込みとして)
かつ `Etale.instOfIsOpenImmersion`(開埋め込みとして)の両方の instance
を持つ——つまり **`C := ∅` は常に有効な `FEt` の余地を作り、`Corr D A B`
は任意の `A B` に対して(空相関により)常に inhabited になってしまう**。
これがあると `Lemma 4.1`(`ZK = D.Ext Z` を要求)は `corrHypInstance4`
上で**文字通り偽**になる(空相関を使えば任意の無関係な `ZK` に対して
反例が作れる)——`Nonempty C` を戻すか、少なくとも `Lemma 4.1` を
証明するときに明示的な追加前提として持ち込む必要がある。★sorry を
埋める前に対処が要る、記録のみ(未修正)。

**発見2(道具の再評価): `.needs` の「rigidity」項目は、おそらく
新しい深い数学(SGA1のπ₁比較定理)ではなく、`AffineTransitionLimit.lean`
に**既にある** Hom-集合の極限安定性の補題群で足りる可能性が高い。**
最初、`.needs` の「étale site の変形に対する剛性」という文言から
SGA1(Grothendieck)の「代数閉体の拡大でπ₁^etが変わらない」という
深い比較定理を疑い、mathlib を検索したが該当なし(`Isom`/`Hom` スキーム
の rigidity は存在しない)。しかし「deformation」という語は文脈上
(EGA IVスタイルの極限による構成)**冪零厚化ではなく有向系の極限**を
指している可能性が高く、その場合に対応するのは
`AlgebraicGeometry.Scheme.exists_hom_hom_comp_eq_comp_of_locallyOfFiniteType`
(`AffineTransitionLimit.lean:686`、2つの段階 `i`・`j` から来る射 `a`・`b`
が極限で一致するなら、ある更に細かい段階 `k` で既に一致する、という
Hom集合の安定性)・`exists_hom_comp_eq_comp_of_locallyOfFiniteType`
(同一段階版)——これは**既に mathlib にあり**、今セッション未使用のまま
残っている道具(前回セッションで「貼り合わせに使う」とだけ記録)。
どちらの読み方が正しいかは実際に組み立てを進めないと確定しないが、
**まず既存の道具で足りると仮定して進め、本当に足りなければそこで
壁を特定する**、という順序が合理的。

**現状の結論**: `Lemma 4.1` の完成には、(a) 発見1の修正
(`Corr`定義への`Nonempty C`の復元、または`Lemma 4.1`への前提追加)、
(b) `c.C` 自体を有限段階 `R` へ貼り合わせながら降ろす構成
(`Scheme.GlueData` を経由、まだ着手していない大きな一手)、
(c) 発見2のHom安定性補題を使った一意性の議論、という3つの残作業がある
——これらは「見つからない数学の壁」ではなく「まだ組み立てていない
エンジニアリング」に分類される点は変わらないが、(a)は原典との
整合性を先に直す必要がある独立した issue として扱う。

### 2026-09-04さらに続報: §1が5/5に到達(`Definition 1.2`を非空虚に実現)

発見1(`Corr`のnonempty脱落)の調査中に見つけた副産物——
`ModularExample.lean` が `.src` を付けずに「正直な限界」として記録して
いた `Definition 1.2`(`Corr.IsTrivial`)の未達成を、`corrHypInstance4`
(`FEt := QcqsFEt`、本物のスキーム射の部分型)で解消した。`corrHypInstance`
(`FEt := PLift ∘ IsFiniteIndexIn`)では証明無関係性により`IsTrivial`が
任意のcorrespondenceについて自動的に真になってしまっていたが、
`QcqsFEt`は本物の射の部分型なので一般には`A≠B`間の`FEt A B`が空にも
複数元にもなり得、`IsTrivial`はcごとに真偽が変わる非自明な述語になる。
`X=Y=C=basePt4`・`α=β=γ=idFEt`という最も単純な例で witness し
(`Found/CorrHyp/TrivialCorrExample.lean`、★sorry無し)、
§1(Definition 1.1–1.5)が**5/5に到達**。全体の集計は7/24→**8/24**。

### 2026-09-04さらに続報: 「命名段」パターンの definition 項目を続けて実現

`Definition 1.2` の解消をきっかけに、§1–§6 の残り項目を`Skeleton`の
statement の**中身**まで見直したところ、`Definition 3.1`(`hyperbolicCore`)
と `Definition 5.2`(`hyperbolicCoreGeneral`)が**どちらも「対象に名前を
与えるだけ」**(非arithmetic性の仮定は定義本体で未使用、と各Skeletonの
docstringが自ら明言している)であることに気づいた。`corrHypInstance`
(`core := id`・`Ext := id`、`Instance.lean` 冒頭で明記済みの placeholder)
で `X := FG_SL2Z`(モジュラー群、非arithmetic性は `MargulisArithmetic`/
`ShimuraArithmetic` が恒偽であることから出る)を使い、両方を実現した
(`HyperbolicCoreExample.lean`・`HyperbolicCoreGeneralExample.lean`、
★どちらも sorry 無し)——`core := id` により witness は定義的に
`FG_SL2Z` 自身に等しいことを、`ModularExample.lean` の透明性のスタイルを
踏襲して docstring に明記した(隠していない)。

一方、`Theorem 6.1`(§6)は`Aut _:=PUnit`・`IsGenericallyScheme _:=True`
という**同じ placeholder パターン**を使えば形式的には閉じるように見えるが、
これは `Definition 1.2` を `corrHypInstance` で claim しなかったのと
**同じ理由で却下**した——`Aut`/`ModuliStack` が本物の自己同型群・
モジュライスタックに結びついていないので、`Theorem 6.1` の主張
(自明な自己同型・correspondence が無い)を実質的に検証したことには
ならない。`Theorem 3.3`(§3)も `Iso X Y := X = Y`(`corrHypInstance`
で本物の等号)により退化の余地が無く、Margulisの理論を要する本物の
有限性の主張として残っている——これらは今回 claim しなかった。

**§1が5/5、§3が2/3、§5が1/7に到達。全体の集計は7/24→10/24。**

### 2026-09-04さらに続報: §2の残り3項目が現状の instance では届かないことを確認

`Definition 1.2`/`3.1`/`5.2` の勢いで §2 の残り(`Proposition 2.4`・
`Theorem 2.5`・`Theorem 2.6`)も同じ手が使えないか検討し、**使えない**
ことを確認した(却下の記録)。

- `Proposition 2.4`(`MargulisArithmetic ↔ ShimuraArithmetic`): `corrHypInstance`
  ではどちらも恒偽なので `False ↔ False` で形式的には閉じるが、これは
  `Definition 1.2` を却下したのと同じ理由(内容を持たない placeholder 同士の
  同値)で claim しない。`corrHypInstance2`(`ShimuraArithmetic` が本物、
  `MargulisArithmetic` は依然恒偽)では `Γ_SL2Z` について
  `False ↔ True` となり**主張が真に偽**——証明できないのは当然で、
  `Definition 2.2`(Margulis-arithmetic の本物の構成)が無い限り届かない。
- `Theorem 2.5`(`Arithmetic ↔ InfinitelyManyCorr`): `Instance.lean` 冒頭の
  docstring がすでに「この instance では成り立たない」と明言済み
  (`InfinitelyManyCorr` 側は非自明になりうるため)——今回の探索で追認しただけ。
- `Theorem 2.6`(有限性): `corrHypInstance` では `Arithmetic` が恒偽なので
  述語 `type X = gr ∧ Arithmetic Γ` が恒偽になり `S := ∅` で**空虚に**閉じて
  しまう——これも `Definition 1.2` 却下と同種の退化なので claim しない。

結論: §2 を実際に前進させるには `Definition 2.2`(代数群 `G` の構成)の
本物の実装が要る——`ShimuraArithmeticData.lean` 自身が「代数群/部分群
スキームの有限性分類が mathlib に無く人年規模」とすでに記録済みの
評価を、今回の探索でも覆せなかった(再確認のみ)。

### 2026-09-04さらに続報: §4の「rigidity」の役割を精緻化

`AffineTransitionLimit.lean` の
`Scheme.exists_hom_hom_comp_eq_comp_of_locallyOfFiniteType` の完全な
シグネチャを確認した:`D : I ⥤ Scheme`・`t : D ⟶ (Functor.const I).obj S`・
`f : X ⟶ S`・`c : Cone D`・`hc : IsLimit c` の下で、`a : D.obj i ⟶ X`・
`b : D.obj j ⟶ X` が「`t` と両立し」かつ「`c` の極限で一致する」なら、
ある共通の細分段階 `k` で両者が一致する、という主張——**これは
「スキーム`C`自体を極限の外から貼り合わせる」道具ではなく、「極限の
頂点上ですでに一致している2つの候補射が、より細かい有限段階で
すでに一致する」という道具**だと判明した。

つまり `Lemma 4.1` の組み立てにおけるこの補題の役回りは、**貼り合わせの
初手(`c.C` 自体を有限段階へ降ろす)ではなく、貼り合わせの整合性検査
(オーバーラップ上で2つの局所モデルが一致することの確認)**である
——初手は`StandardEtalePair`の道具(`exists_finite_standardEtaleCover`+
`standardEtalePairRingBaseChange`)がすでに担っている。これで
「rigidity は既存道具で足りるか」という前回の疑問に、**役割としては
Yes(ただし貼り合わせ段でのみ使う)**という、より正確な答えが出た。

**次の具体的な一手**(まだ未着手): `exists_extDiagram_finite_affine_descent`
が返す1つのアフィン片について、`c.α` をその片へ制限した先の
`Algebra.Etale` 環(`Etale.algebraEtale_appLE`)を
`exists_finite_standardEtaleCover` で標準エタール片に細分し、各片を
`standardEtalePairRingBaseChange` で有限段階まで降ろす——この「1片の
降下」を実際に組み立てるには、`Γ(その片, ...)` が `Γ(有限段階の対応する片,
...)` の `K` への base change であることを、`isLimit_extCone` の
Γ(大域切断)版ではなく**任意のアフィン開集合版**として持ち上げる必要が
あり、まだ書いていない——次のセッションの最初の一手として記録する。

### 2026-09-04さらに続報: 「1片のΓをbase changeとして書く」ための鍵の道具を特定

上記の欠けている環同型(`Γ(U_j,U_j) ≅ Γ(V_j,V_j) ⊗_R K`、`U_j` は
`(extCone X).π.app R` による `V_j` の preimage)を実際に作る道具として
**`AlgebraicGeometry.pullbackSpecIso`**(mathlib 既存)を特定した:
`pullback (Spec.map (algebraMap R S)) (Spec.map (algebraMap R T)) ≅
Spec (S ⊗[R] T)`——アフィン同士の pullback がテンソル積の Spec そのもの
であることを直接与える。`AffineTransitionLimit.lean` 側は「Hom集合・
大域切断の安定性」しか持たないので、この橋渡し自体は**mathlibのpullback
理論から**(`AffineTransitionLimit.lean`からではなく)持ってくる必要が
あった、という発見。

残る作業: `U_j = π.app R ⁻¹ᵁ V_j` が実際に `V_j ×_{Spec R} Spec K`
(pullback square)であることを`extConePi`/`extDiagram`の構成
(`Limits.pullback.map` ベース)から示し、`V_j` がアフィン
(`hVprop j`)なので `pullbackSpecIso` を適用して `Γ(U_j,U_j) ≅
Γ(V_j,V_j) ⊗_R K` を得る——`extConePi_app_fst`/`_snd`(`ExtLimit.lean`
既存)がこの pullback square の脚をすでに与えているはずなので、
そこから `IsPullback` の四角を組み立てる作業になる。まだ書いていない。

### 2026-09-04さらに続報: 「`Ext X` の射影 = 有限段階への base change」を完成 ★★★★★

上記の残作業のうち、**pullback square の存在そのもの**を完成させた
(`ExtLimit.lean`、`extConePi_app_eq`、★sorry無し・標準3公理のみ)。

鍵は mathlib の `CategoryTheory.Limits.pullbackLeftPullbackSndIso`
(`pullback (pullback.snd f g) g' ≅ pullback f (g' ≫ g)`、pullback の
pasting)——`f := X.hom`・`g := (toSchemeDiagramOver.obj R).hom`・
`g' := phiR R`(`Spec K → Spec R` の遷移射、`specKCone` の射影)とおき、
`phiR R ≫ (toSchemeDiagramOver.obj R).hom = toBaseK`(`phiR_comp`、
`specKConeOver` の cone 条件そのもの、`isLimit_specKConeOver` の証明が
特定の代表元 `R0` だけで使っていた事実の一般形)で右辺を `toBaseK`
に付け替えると(`pullback.congrHom`)、**`Ext X` の台 = 有限段階 `R`
のファイバー積 `P_R` を `Spec K → Spec R` に沿って base change したもの**
という同型 `extConeIso` が得られる。さらに**`extConePi X .app R`
(`Ext X` の `R` への射影)が、この同型の下で「`P_R` への射影」そのもの**
であること(`extConePi_app_eq`)を、`extConePi_app_fst`/`_snd`(定義)と
`pullback.hom_ext` の照合で示した。

配管の教訓(`tools/lean-idioms.md` 第25項に追記): 第22項の「`have h' :
<明示型> := h` で defeq を迂回する」技を **`Type` の値(`Hom`)**に適用する
ときは `have` ではなく **`set`(型注釈つき)** を使う必要がある——`have`
は`Prop`では問題にならないが、`Type`の値の中身を後続の項から見えなく
してしまうため。

**これで「1アフィン片の降下」に必要な純粋に幾何学的な部分が完成した**
——残るのは `pullbackSpecIso`(mathlib既存)を`V_j`(アフィン)に適用して
`Γ(U_j,U_j) ≅ Γ(V_j,V_j) ⊗[R] K` という**環レベルの結論**を引き出す
最後の一手(`extConeIso`の`V_j`への制限、まだ未着手)。

### 2026-09-04さらに続報: 残る環化の手順を特定(未着手)

`pullbackSpecIso`を`V_j`(`(extDiagram X).obj R`のアフィン開)に実際に
適用するための残り手順を特定した:
1. `Q.fst ⁻¹ᵁ V_j`(`Q := pullback (pullback.snd X.hom (toSchemeDiagramOver.
   obj R).hom) (phiR R)`)を`pullbackRestrictIsoRestrict`(mathlib既存)
   + `pullbackRightPullbackFstIso`(pullback の pasting、既出)で
   `pullback (V_j.ι ≫ pullback.snd X.hom (toSchemeDiagramOver.obj R).hom)
   (phiR R)` と同一視する(pullbackSymmetryが要るかもしれない)。
2. `V_j.ι ≫ pullback.snd X.hom (toSchemeDiagramOver.obj R).hom : V_j ⟶
   Spec R` を、`V_j`がアフィン(`hVprop j`)なので**`Spec.map`の形に書く**
   ——鍵は`AlgebraicGeometry.arrowIsoSpecΓOfIsAffine`(mathlib、アフィン
   同士の射は`Arrow.mk f ≅ Arrow.mk (Spec.map f.appTop)`)+
   `IsAffineOpen.isoSpec`(`V_j ≅ Spec Γ(V_j,V_j)`)。
3. `pullbackSpecIso`を適用して`Γ(U_j,U_j) ≅ Γ(V_j,V_j) ⊗[R] K`を得る。

まだ未着手(道具の在庫確認のみ)。次のセッションの最初の一手として記録。

### 2026-09-04さらに続報: 手順2を実測、かつ残作業の規模がもう1段階あることが判明

手順2(`V.ι ≫ pullback.snd ... : V ⟶ Spec R` を `Spec.map` の形に書く)を
実際に試し、`arrowIsoSpecΓOfIsAffine`(`g : X ⟶ Y` に `[IsAffine X]
[IsAffine Y]` があれば `Arrow.mk g ≅ Arrow.mk (Spec.map g.appTop)`)が
`V`(`hV : IsAffineOpen V` から `IsAffine ↥V`)と`(toSchemeDiagramOver.obj
R).left`(新たに`toSchemeDiagramOver_obj_isAffine`で`IsAffine`を確認、
`ExtLimit.lean`に追加・★sorry無し)の両方に対して直接使えることを確認
——`Arrow.mk`の同型から`.hom.left`/`.hom.right`/`.hom.w`で実際の同型
射・可換四角を取り出せることも確認した。手順1-3の**道具はすべて実在し
使える**ことが実測で裏付けられた。

一方、この過程で**残作業がもう1段階あること**に気づいた——
`Γ(U_j,U_j) ≅ Γ(V_j,V_j) ⊗[R] K` を得たあとも、`exists_finite_
standardEtaleCover`で得られる「`Γ(U_j,U_j)` 上の標準エタール対」を
**この局所環 `Γ(V_j,V_j)` 上まで降ろす**という、これまでの
`FieldLimit.lean`(`exists_fg_subalgebra_standardEtalePair`、`k=ℚ`・
`K=ℝ`という原始の対でのみ構築済み)を**局所環 `Γ(V_j,V_j)` を新しい
"底"とする版へ一般化する**作業がなお必要——手法(多項式係数の有限生成
部分環への還元)自体は使い回せる見込みだが、宣言としては未着手。
「見つからない数学の壁」ではない点は変わらないが、残工程が当初の想定
より1段階多いことを正直に記録する。

### 2026-09-04さらに続報: 訂正——上記の「一般化」は実は不要、既に汎用だった

上記の懸念を検証しようと `exists_fg_subalgebra_standardEtalePair_map`
(`FieldLimit.lean`)の実際の型を読み直したところ、**`{k K : Type*}
[CommRing k] [CommRing K] [Algebra k K]` という完全に一般の型で
すでに書かれており、`k=ℚ`・`K=ℝ` に限定されていなかった**——先の
「局所環版への一般化が要る」という評価は誤りで、`k := Γ(V_j,V_j)`・
`K := Γ(U_j,U_j)` を直接代入するだけで適用できる。この点は訂正する
(過大評価を修正、朗報)。

ただし、代わりに**別の、より本質的な難点**が見えてきた:
`exists_fg_subalgebra_standardEtalePair_map` が返す `R'`(標準エタール対
が実際に住む有限生成 `Γ(V_j,V_j)`-部分環、`K=Γ(U_j,U_j)=Γ(V_j,V_j)⊗_R K`
の中)は、テンソル積の中の抽象的な部分環であり、**それ自体が「`Spec R''`
という、`FgSubalgebra ℚ ℝ` の中のさらなる細分段階」に自動的には対応
しない**——対応させるには、`R'` を生成する有限個の元(テンソルの和)を
さらに `ℚ`-部分環の言葉に翻訳し直す、別の「有限生成部分環への還元」の
議論が要る。さらに、**複数のアフィン片 `V_j` それぞれで得られる細分段階
`R'_j` を、全体で共通の1つの細分段階 `R''` へ合流させる**という、
`AffineTransitionLimit.lean` のHom安定性(`exists_hom_hom_comp_eq_comp_
of_locallyOfFiniteType` 等、前々回セッションで役割を特定済み)を実際に
呼び出す段が必要——ここが「rigidity」の出番だと再確認できた。

**現状のまとめ**: 個々の道具(環化・標準エタール対の降下・Hom安定性)は
すべて特定済みで実在するが、**それらを1つのアフィン片の中で組み合わせ、
かつ複数の片を横断して1つの共通段階へ合流させる**という、2段階の
組み立て作業が残っている。数学的な内容としては段々と輪郭がはっきりして
きているが、Lean宣言としては依然未着手。

### 2026-09-04さらに続報: 「合流」段に潜む本物の数学的難点(flatness)を発見

上記「複数の片を横断した合流」を実際に組み立てようとして、**見た目の
配管ではなく本物の数学的な難点**を発見した——`Γ(V_j,V_j) ⊗[R] ℝ` の
中の有限個の元(標準エタール対の係数)を、有限生成部分環 `R'' ⊇ R` を
使って `Γ(V_j,V_j) ⊗[R] R''` まで引き戻すこと自体は易しい(各元を
純テンソルの有限和として書き、現れる実数係数を全部 `R` に添加すれば
`R''` が作れる、既存の `exists_fg_subalgebra_polynomial_family` と同種の
議論)。**しかし**、引き戻した先で `StandardEtalePair` の定義方程式
(`f' * p₁ + f * p₂ = g^n`)が**実際に成り立つこと**を示すには、
`Γ(V_j,V_j) ⊗[R] R'' → Γ(V_j,V_j) ⊗[R] ℝ` が**単射**である必要がある
——これは元の`exists_fg_subalgebra_standardEtaleCond`が
`Polynomial R`(`R`上**自由**加群なので `algebraMap` の単射性から
`Polynomial.map`の単射性が自動的に従う)を使っていたのと違い、
`Γ(V_j,V_j)` は一般には `R` 上**自由でも平坦でもない**(スキーム `X`
の構造射が平坦とは限らないため)——**このままでは単射性の根拠が無い**。

標準的な解決策は EGA IV の**generic flatness**(有限型加群は、底環の
ある稠密開集合上では平坦になる)を使って `R` をさらに細分し、その先で
`Γ(V_j,V_j)` が平坦になるようにしてから議論する、というもの——
これは「配管」ではなく**もう1つ本物の代数幾何の定理**であり、
`.cache/mathlib-index.txt` を `generic.*flat`/`genericFreeness` で
検索したが**該当無し**——2026-09-04時点でmathlibに generic flatness
は無いと見られる(次セッションで exact? 等による直接確認が要る)。
これで `Lemma 4.1` の残作業には、(a) 合流の配管、(b) generic
flatness(または平坦性を回避する別の議論)、という**もう1段掘り下げた
本物の数学**が要ることが分かった——★依然「壁ではなく道」だが、
過去に見積もった「2段階」よりもう少し深い。

### 2026-09-04さらに続報: (b)にはさらに前提が要ることが判明——`Space`の有限型性

generic flatness(EGA IV 6.9.1、Eisenbud定理14.4 "generic freeness")
自体の標準的な前提を確認した——`R` はネター整域(またはreduced)、
かつ加群/多元環は `R` 上**有限表示**(finitely presented、少なくとも
finitely generated as an algebra)である必要がある。ところが
`corrHypInstance4` の `Space := QcqsSpace` は **qcqs(コンパクト・
準分離)しか要求しておらず、有限型(locally of finite type)を要求して
いない**——一般の `X : QcqsSpace` に対して `Γ(V_j,V_j)` が `R` 上
有限生成である保証がなく、generic flatness が(たとえ形式化できても)
**そのままでは適用できない**。

これは新たな発見であり、`Lemma 4.1` を完全に誠実に証明するには
**`Space` をさらに「qcqs かつ有限型」に絞った新しい instance**
(`QcqsSpace` の上にもう1段、`corrHypInstance5` 相当)が要ることを
意味する——双曲曲線は本来有限型なので原文の意図には忠実(§4の
`corrHypInstance3→4` の pivot と同種の、原文に忠実な絞り込み)。

**現状の総括**: `Lemma 4.1` の完全な証明には、(1) `Space` を有限型
qcqsへさらに絞る、(2) generic flatness(またはその代替)を形式化する、
(3) 複数アフィン片の細分段階を合流させる配管、(4) `GlueData` での
貼り合わせ、という**4段階の残作業**があることが、今回の掘り下げで
明確になった。個々の段は(1)は既存パターンの反復、(4)は在庫確認済み
だが、(2)は独立した本物の代数の定理として新規に要る——これは
1回のセッションで詰め切れる規模を超えており、複数セッションにまたがる
取り組みとして継続する。★★引き続き「壁ではなく道」——今回のセッション
で「道の残り距離」が具体的に測れるようになったこと自体が前進である。

### 2026-09-04さらに続報: 大きな前進——generic flatnessを完全に回避する戦略が実現 ★★★★★

(2)(generic flatness)を実際に回避する**戦略の転換**を見つけ、
`ExtLimit.lean` に `piecePullbackIso` として完成させた(★sorry無し・
標準3公理のみ、コミット `80eb752a`)。

**転換の核心**: これまでの戦略(`exists_extDiagram_finite_affine_
descent`)は、`P_R`(有限段階 `R` のファイバー積 `X.left ×_k Spec R`)の
**任意の**アフィン開`V_j`から出発していた——この`V_j`の切断環`Γ(V_j,V_j)`
は`R`上一般には平坦とは限らず、これが generic flatness を要求する
根本原因だった。代わりに、**`X.left` 自身の(`R` に依らない)アフィン
開被覆 `{U_i}`** をそのまま `Ext X = X.left ×_k Spec K` へ base change
すれば、各片は必ず `Spec(Γ(U_i,U_i) ⊗[ℚ] ℝ)` になる——**`ℚ` は体**
なので `Γ(U_i,U_i) ⊗[ℚ] -` は**自動的に完全関手(平坦)**になり、
generic flatness が構造的に不要になる。

**証明の骨格**(`piecePullbackIso`): `pullbackRestrictIsoRestrict`+
`pullbackSymmetry`+`pullbackRightPullbackFstIso`(pullback の pasting、
`extConeIso`で確立した技法の再利用)で `(pullback.fst X.hom toBaseK)
⁻¹ᵁ U` を `pullback (U.ι ≫ X.hom) toBaseK` と同一視し、`Spec.preimage`
(`Spec.homEquiv`の逆写像)で作った環準同型`pieceRingHom`(定義方程式
`pieceRingHom_spec`)+`pullbackHomIsoLeft`(mathlibに無かった「片方の
脚を同型で付け替えても pullback は変わらない」、`pullback.map_isIso`
から補った)で `Spec.map` 形に付け替えてから `pullbackSpecIso` を適用。

**これで残作業4段階のうち(2)generic flatnessが不要になった**——
残るのは(1)`Space`の有限型性(この`piecePullbackIso`アプローチでは
`Γ(U_i,U_i)`が`ℚ`上有限生成であることまでは要らない、`X.left`のみが
compact/qcqsであれば足りる可能性が高い、要再検証)・(3)複数`U_i`の
細分段階の合流(こちらは`ℚ`上のテンソル分解なので既存の
`exists_fg_subalgebra_polynomial_family`パターンがそのまま使える見込み)
・(4)`GlueData`貼り合わせ、の3段階に縮小した。★★これは今セッション
の中でも特に価値の高い発見——「壁だと思っていた場所に、実は回り道が
あった」という好例。

### 2026-09-04さらに続報: (3)「複数片の合流」に要る単射性も確認済みと判明

(3)(複数`U_i`の細分段階の合流)の核心は「`Γ(U_i,U_i) ⊗[ℚ] R'' →
Γ(U_i,U_i) ⊗[ℚ] ℝ`(`R'' ⊆ ℝ` が細分段階)が単射である」ことだが、
これは**`ℚ` が体なので `Module.Flat ℚ A` が任意の `ℚ`-代数 `A` で
`infer_instance` 一発で通る**(体上の加群は常に平坦)ことと、
`Module.Flat.iff_rTensor_preserves_injective_linearMap`(mathlib既存、
平坦性 ↔ テンソルが単射性を保つ)を組み合わせるだけで**既に確認できた**
(小さな例で実測、`ExtLimit.lean`にはまだ書いていない)——このため
(3)は本当に「既存パターンの反復」で閉じる見込みが高いことが裏付けられた。

**§4の現状**: `Lemma 4.1`の完全な証明に要る道具は、(2)generic flatness
が完全に不要になったことと、(3)の単射性が確認できたことにより、
**「見つからない数学」は無くなった**——残るのは(1)`Space`の有限型性の
検証、(3)の実際の組み立て(有限個の元のテンソル分解→係数の収集→
`exists_fg_subalgebra_mem_finset`適用、という具体的な配管)、
(4)`GlueData`貼り合わせ、という**エンジニアリングのみ**。数学的な
見通しとしてはこのセッションで大きく前進した。

### 2026-09-04さらに続報: (3)を実際に完成——「1アフィン片の降下」の核が揃った ★★★★★

上記(3)を実際に`FieldLimit.lean`へ書き切った。`exists_fg_subalgebra_
tensor_finset`(`A ⊗[ℚ] ℝ` の有限個の元は、ある `R : FgSubalgebra ℚ ℝ`
を使って `A ⊗[ℚ] R.1` からの像として同時に書ける)→
`exists_fg_subalgebra_tensor_polynomial_family`(多項式の係数版)→
`tensor_map_injective_of_flat`(`Module.Flat ℚ A` が `infer_instance`
一発で通ることから単射性を得る)→`exists_fg_subalgebra_tensor_
standardEtaleCond`→**`exists_fg_subalgebra_tensor_standardEtalePair`**
(`StandardEtalePair (A ⊗[ℚ] ℝ)` は、ある `R` を使って
`StandardEtalePair (A ⊗[ℚ] R.1)` の像として書ける)という、既存の
`exists_fg_subalgebra_standardEtalePair`(一般の`k K`)と**完全に並行な
構成**を、`k := ℚ`(体)・`K := A ⊗[ℚ] ℝ` という特殊化で書き切った。
★**すべてsorry無し・標準3公理のみ**。

**これで「`piecePullbackIso`(`Ext X` のアフィン片は `Spec(Γ(U,U)⊗[ℚ]ℝ)`)
+ この一連の降下補題」を組み合わせれば、`Lemma 4.1`の「1アフィン片の
降下」に要る核となる部品がすべて揃った**——`c.α`(相関の有限エタール射)
をアフィン片へ制限 → `Etale.algebraEtale_appLE` で `Algebra.Etale` →
`exists_finite_standardEtaleCover` で標準エタール片に細分 →
`exists_fg_subalgebra_tensor_standardEtalePair` で有限段階 `R` へ降下、
という**一直線の道筋**が(数学的には)完成した。残るのは、この道筋を
実際に1つの `Found/` 宣言として組み立てる作業と、複数のアフィン片を
横断した整合性(rigidity)・`GlueData` 貼り合わせ。

### 2026-09-04さらに続報: base change 復元まで完成——1片降下の核心が完全に揃った ★★★★★

`exists_fg_subalgebra_tensor_standardEtalePair`(降下)と
`standardEtalePairRingBaseChange`(base change 可換性、既存)を組み合わせ、
**`exists_fg_subalgebra_tensor_standardEtalePair_baseChange`**を完成
(`FieldLimit.lean`、★sorry無し・標準3公理のみ)。「`StandardEtalePair
(A⊗[ℚ]ℝ)` は、ある有限段階 `R` 上の `StandardEtalePair` の base change
として**文字通り一致する**」(`P₀.Ring` を `A⊗[ℚ]ℝ` へ係数拡大すると
`P.Ring` に同型)ことまで保証した——`P₀.map(algebraMap...) = P` という
**構造体そのものの一致**(`f`・`g` の一致から `monic_f`・`cond` は Prop
無関係性で自動的に従う、`cases`+`subst`で示した)を経由。

**これで`Lemma 4.1`の「1アフィン片の降下」に要る核心が完全に揃った**:
`piecePullbackIso`(`Ext X`のアフィン片`=Spec(Γ(U,U)⊗[ℚ]ℝ)`)+
`Etale.algebraEtale_appLE`(スキーム→環の橋渡し)+
`exists_finite_standardEtaleCover`(étale-locus細分)+この補題(有限段階
への降下+base change復元)、という一直線の道筋が**数学的にもLean宣言
としても揃った**。残るのは(a)これらを実際に1つの`Found/`宣言として
組み立てる作業、(b)複数のアフィン片(有限個のstandard-étale片、
有限個の`U_i`)を横断した細分段階の合流、(c)`GlueData`貼り合わせ、
という3段階のエンジニアリングのみ。

### 2026-09-04さらに続報: (b)「合流」に要る道具も出揃った

`exists_fgSubalgebra_upperBound`(有限個の`FgSubalgebra k K`は共通の
上界を持つ、`IsDirected`をFinset上の帰納法で拡張)・
`StandardEtalePair.map_map`(`(P.map φ).map ψ = P.map(ψ.comp φ)`、
`f`・`g`の一致から構造体の一致を出す)を`FieldLimit.lean`に追加
(★どちらもsorry無し)。これで(b)「ある段階`R_i`で得た局所モデルを、
より粗い共通段階`R'`(有限個の`R_i`の上界)へ移送してもbase changeが
変わらない」ことを示す道具が揃った——`exists_fg_subalgebra_tensor_
standardEtalePair_baseChange`の出力(段階`R_i`ごとに異なる)を、
`exists_fgSubalgebra_upperBound`で得た共通`R'`へ`StandardEtalePair.
map_map`で移送すれば、有限個の標準エタール片をすべて**同じ`R'`**の
上で扱えるようになる。

**§4の現状**: (a)実際の組み立て・(c)`GlueData`貼り合わせを除き、
`Lemma 4.1`の証明に要る個々の数学的補題は**すべて完成した**。残るのは
これらを実際に1本の証明として結線するエンジニアリング(かつ`Corr`の
nonempty脱落・`Space`の有限型性という、独立に対処が要る2つの前提の
整理)のみ。

### 2026-09-04さらに続報: スキーム→環の橋渡しを完全に結線 ★★★★★

(a)(実際の組み立て)を一段進め、`c.α : C ⟶ Ext X`(有限エタール)を
`X.left` 由来のアフィン片へ制限すると `exists_finite_standardEtaleCover`
・`exists_fg_subalgebra_tensor_standardEtalePair_baseChange` が**直接
読める形**(`Γ(U,U)⊗[ℚ]ℝ` 上の `Algebra.Etale`)になる、というところまで
`ExtLimit.lean` に結線した(`piece_algebraEtale_tensor`、★sorry無し)。

鍵は3つの補題の合成: `piece_algebraEtale`(スキーム→環、環は
`Γ((ExtF.obj X).left,V)`のまま)→`pieceRingEquiv`(`piecePullbackIso`
を`Scheme.Opens.topIso`・`Scheme.Γ.mapIso`・`Scheme.ΓSpecIso`で繋いで
環の同型に変換)→`algebraEtale_transport`(`RingHom.Etale.respectsIso`
で底環の同型に沿って`Algebra.Etale`を輸送、mathlib既存の道具の新しい
組み合わせ)。

**これで「`c.α` → アフィン片への制限 → `Algebra.Etale (Γ(U,U)⊗ℚℝ) ...`
→ `exists_finite_standardEtaleCover` → `exists_fg_subalgebra_tensor_
standardEtalePair_baseChange` → 有限段階への降下」という、`Lemma 4.1`
の「1アフィン片の降下」の**全行程が実際に1つのLean宣言の連鎖として
結線可能な状態**になった。残るのは、この連鎖を実際に呼び出して
1つの完成した宣言にまとめる作業(標準エタール被覆の有限個の片を横断
した合流を含む)と、複数の`U_i`・`GlueData`貼り合わせ。

### 2026-09-04さらに続報: 「合流」の最後の道具(移送)が完成 ★★★★★

`exists_fg_subalgebra_tensor_standardEtalePair_promote`を完成
(`FieldLimit.lean`、★sorry無し)——有限段階`R`上の降下`P₀`を、より粗い
共通段階`R'`(`R≤R'`)へ移送しても base change が変わらないことを保証。
`StandardEtalePair.map_map`+`Algebra.TensorProduct.map_comp`(mathlib、
`⊗`の関手性)を組み合わせるだけ。

**これで「1つのアフィン片`U`の上で、標準エタール被覆の有限個の片
(それぞれ異なる降下先`R_i`を持つ)を横断して1つの共通段階`R'`へ合流
させる」ために要る道具が完全に揃った**——各`i`について`exists_fg_
subalgebra_tensor_standardEtalePair_baseChange`で`R_i`・`P₀_i`を得て、
`exists_fgSubalgebra_upperBound`で共通の`R'`を得て、`exists_fg_
subalgebra_tensor_standardEtalePair_promote`で各`P₀_i`を`R'`へ移送する
——という3段の組み合わせで、1つのアフィン片`U`の中の**すべての**標準
エタール片を同じ有限段階`R'`の上で扱えるようになる。

**§4の現状(総括)**: `Lemma 4.1`の証明に要る数学的な補題・道具は
(1)`Space`の有限型性の検証という前提整理を除き**すべて完成した**。
残る作業は純粋なエンジニアリング(実際に結線して1つの宣言にまとめる、
複数の`U_i`を横断した合流、`GlueData`貼り合わせ、`Corr`のnonempty
脱落の手当て)のみであり、必要な数学はすべて手元に揃っている。
