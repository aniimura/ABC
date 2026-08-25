# 配管の対策集 —— この codebase で繰り返し出る失敗形と、効いた直し方

★これは**数学の難しさではなく、エラボレータとの戦い方**の記録である。
同じ穴に何度も落ちたので、1 箇所にまとめる。

★★**運用**: 新しい失敗形に当たったら、直すと同時にここへ 1 行足すこと。
「前にも見た」と思ったらまずここを引く。

---

## 1. `instances 透明度で型が合わない` —— 最頻出

**症状**: `rw` / `simp` / インスタンス探索が
「`Did not find an occurrence of the pattern ...`」と言って落ちる。
目標にはその式が**そのまま見えている**。エラー末尾に

```
Note: The target expression is not type-correct under the `instances` transparency level
```

が付く。原因は `Under.right W` / `Z.unop.left.obj` / `WideSubcategory` /
`PfRootObj` の射影が、`instances` 透明度では展開されず型が合わなくなること。

**効いた直し方(上から順に試す)**

| 手 | 例 |
|---|---|
| ★`rw` をやめて **`Eq.trans` / `congrArg` の項**で繋ぐ | `(Category.assoc _ _ _).trans (h.trans (Category.assoc _ _ _).symm)` |
| ★`congrArg` には**関数の型を明示**する | `congrArg (fun t : X ⟶ Y => t ≫ f) h` |
| ★`show` で**きれいな型**に言い換えてから触る | `show HomBirat.mk (biratPfIdx …) (rootMap … ≫ …) = _` |
| 射影を**名前つき `def`** に逃がす | `rtObjPf` / `pfDown`(既存の実例) |
| 添字対象を `idxBiratMk P G a hc hs` の形で**受け取る** | 型が素の `C` の対象になる |

## 2. `cancel_epi` / `cancel_mono` がインスタンスを見つけない

**症状**: `have hep : Epi f := …` を置いたのに
`(cancel_epi _).mp h` が `failed to synthesize Epi f` で落ちる
(`haveI` にしても同じ)。

**直し方**: **構造体の射影を直接使う**。

```lean
have hep : Epi f := P.totEpiC _ _ f
exact hep.left_cancellation _ _ h        -- cancel_epi の代わり
exact hmono.right_cancellation _ _ h     -- cancel_mono の代わり
```

★`Epi` / `Mono` は `left_cancellation` / `right_cancellation` を
フィールドに持つので、インスタンス探索を経由せずに済む。

## 3. `ℕ≥1` が**依存位置**(対象の中)に現れる

**症状**: `⟨A, n⟩ : PfRootObj P F` のように次数が対象の一部になっていると、
「次数は命題として等しいが項としては違う」2 つの構成が**同じ項にならない**。
`rw` は `motive is not type correct` で落ちる。

**直し方**: ★★**次数を仮引数 `K` に出し、等式を仮引数で受け取る**。

```lean
def fooK (W : Idx …) (K : ℕ+) (hK : degFr W = K) : … := …   -- 一般形
theorem fooK_eq … : fooK W K hK … = foo W … := by subst hK; rfl
```

`K` は**変数**なので `subst` が効き、証明部分は `Prop` なので
**証明が違っても同じ項**になる。
実例: `biratPfIsoA'` / `biratPfMk'` / `biratPfMk'_eq`(`Prop55BiratPf.lean`)。

★**規約**: 依存位置に `ℕ≥1` を持つ定義を書いたら、
**その場で primed 版(次数を仮引数に出した版)も並置する**。
後から足すと呼び出し側を全部書き直すことになる。

## 4. `subst` を使うために、等式の片側を**変数**にする

**症状**: `h : f x = c.field` の `c.field` は変数ではないので `subst` できない。

**直し方**: 補題の仮引数に `δ` を取り、`hδ : f x = δ` の形で受ける。
呼び出し側で `δ := c.field` を渡せばよい。
実例: `biratPfHom_surj_mk`(`δ` と `φ` を仮引数に出してから `subst`)。

★構造体の **eta** が効くので `idxBiratMk P G T.unop.hom.hom _ _ = T` は **`rfl`** である
(`Over` / `WideSubcategory` / `Discrete PUnit` / `Opposite` すべて eta を持つ)。
これで「一般の添字対象」を「構成子の形」に**無料で**言い換えられる。

## 5. `rw` の末尾 `rfl` が閉じない

**症状**: 目標が見た目 `X * Y = X * Y` なのに `unsolved goals`。

**直し方**: `rw` の自動 `rfl` は `reducible` 透明度なので、
`exact rfl` / `exact mul_comm _ _` / `exact congrArg (· * y) h` と**項で閉じる**。

## 6. `𝟙` がどの圏の恒等射か決まらない

**症状**: `PfCat P F` は定義上 `C` そのものなので、
`𝟙 (X : PfCat P F)` と書いてもインスタンスは `C` の方が選ばれる。
型注釈を対象に付けても直らない。

**直し方**: **恒等射を像として書く** —— `toHomPf (𝟙 X)` / `toHomBirat (𝟙 X)`。
関手の `map_id` が `rfl` なので、これは本当に恒等射である。
実例: `biratPfHom_id`。

## 7. 巨大な構造体を型の中に直接書かない

**症状**: `Frobenioid (pfRootPre P F)` を型の中に書くと
`whnf` が 200000 heartbeats を超える(**実測 73 秒で timeout**)。

**直し方**: **仮引数で受ける**(`(Gpf : Frobenioid (pfRootPre P F))`)。0.2 秒になる。

## 8. `Under.isoMk` / 無名構成子が未簡約の型に当たらない

**直し方**: 同型を**型注釈つきの `let`** で別に組んでから使う。

```lean
let e : W' ≅ d.right := { hom := …, inv := …, … }
```

## 9. 宇宙変数の使い回し

**症状**: 始域と終域で別の圏を扱う補題で `u2 v2` を使い回すと当たらない。

**直し方**: `universe uu1 vv1 uu2 vv2` と**新しく宣言する**。
実例: `iso_unique_of_rigid`(`Cor54Rigid.lean`)。

## 10. ファイル書き込みの事故

* Lean のファイル内容は **`Write` ツールで書く**。
  python の heredoc に `𝒞` のようなサロゲート対のエスケープを入れると
  Windows の既定エンコーディングで `UnicodeEncodeError` になり、
  **ファイルが途中まで書かれて壊れる**(実際に 772 行失った)。
* `git commit -m` にバッククォートを入れない。`git commit -F -` ＋ heredoc を使う。
* `lake build` は `lean/` から、`node tools/check.mjs` は**リポジトリ直下**から。

## 11. doc コメントと宣言の間に `set_option ... in` を置けない

**症状**: `unexpected token 'set_option'; expected 'lemma'`。
MCP(`lean_check`)では通るのに `lake build` で落ちるのでたちが悪い
——MCP へ投げる断片には doc コメントを付けないことが多いため。

```lean
/-- ★説明 -/
set_option maxHeartbeats 1000000 in
theorem foo : … := …          -- ✗ パースエラー
```

**直し方**: `set_option ... in` を **doc コメントより前**に書く。

```lean
set_option maxHeartbeats 1000000 in
/-- ★説明 -/
theorem foo : … := …          -- ✓
```

`variable (P F) in` / `omit … in` も同じ。★**属性 `@[simp]` は逆で、doc の後**である。

## 12. 「包み `def` の射影」が `rw` の照合を止める

**症状**: 補題は `X.obj` について述べているのに、目標には
`(scaleRootObj k X).obj` と出ていて `rw` が当たらない。定義上は同じ
(`scaleRootObj k X := ⟨X.obj, k * X.root⟩`)だが**構文が違う**。
1 と同じ `instances` 透明度の注記が付くことが多い。

**直し方**: 触る直前に**包みだけ `unfold`** する。

```lean
rw [hR, hL]
unfold scaleRootHom scaleRootObj   -- ← これで (scaleRootObj k X).obj が X.obj になる
rw [hf, hg]
```

★★`simp only [scaleRootObj]` でも同じだが、`unfold` の方が
「包みを 1 枚剥いだだけ」であることが読み手に分かる。
★依存位置の `ℕ≥1`(3 番)と合わせて、**包みの中の根は必ず項として書き下せる形にしておく**。

## 13. `IsIso` のインスタンスが `haveI` で登録されない

**症状**: 仮定に `ha : IsIso (P.Base a)` があり `haveI := ha` もしたのに、
`IsIso.hom_inv_id` などが `failed to synthesize instance IsIso (P.Base a)` と言う。
★名前をつけて `haveI hA : IsIso (P.Base a) := ha` としても直らないことがある
(`rw` で式を書き換えた**後**の目標では、射の項が `instances` 透明度で
別物に見えているため)。

**直し方**: **インスタンスを明示的に渡す**。

```lean
exact @IsIso.hom_inv_id _ _ _ _ (P.Base a) ha       -- ✓ 確実
-- exact IsIso.hom_inv_id _                          -- ✗ 合成に失敗しうる
```

★`inv` を含む等式を述べるときも `@inv _ _ _ _ f ha` と書いておくと、
呼び出し側でインスタンスがずれない。`IsIso` は `Prop` なので
**証明無関係により、どの証明を渡しても同じ項**になる。

---

## 検査器のキャッシュ(2026-08-21)

`node tools/check.mjs` の時間はほぼ全部が `pdftotext` の呼び出しだった。
`.cache/pdf-pages.json` に**ページ本文を跨いでキャッシュ**するようにした
(**55 秒 → 7 秒**)。

★鍵は `PDF のパス # ページ # mtime # size` に加えて
**`check.mjs` 自身のハッシュ**を含む —— `normalize` / `squash` / `PDF_MODES` を
触ったら必ず外れる。これを忘れると
「正規化を変えたのに古いテキストで通る」という**器具の穴**になる。

---

## 14. 型の同義語(`PfCat P F := C`)を跨ぐと `rw` が当たらない(2026-08-24)

**失敗形**: `PfCat P F` は `C` の型の同義語で、`pfDiv` などは
暗黙引数を `{A B : PfCat P F}` で取る。一方 `toHomPf ψ` は `ψ : A ⟶ B`
(`A B : C`)から作るので、**同じ項でも暗黙引数の書かれ方が違う**。

```
-- 目標側:   pfDiv (A := pfDown P F ((toPfCat P F).obj X)) …
-- 補題側:   pfDiv (A := rtObj P F (pfObjDown P F X) 1) …
```

`rw` は `instances` 透明度でしか合わせないので、
**defeq なのに「パターンが見つからない」**になる。

★`(toPfCat P F).obj A = A` の simp 補題を足しても、今度は逆向きに
「`C` の項が `PfCat P F` を期待されている」でずれる。

**対処(測定済み)**:
* 同義語の**両向き**に名前付きの `abbrev` を置く(`pfObjDown` / `pfUp`)。
* それでも駄目なら `rw` をやめて `Eq.trans` / `congrArg` の**項の側**で組む
  (`exact` は既定透明度で合わせるので通る)。
* 圏の合成 `≫` を `compPf` に開く橋 `pfCat_comp_eq : f ≫ g = compPf P F f g := rfl`
  を用意しておくと `pfDiv_comp` 等が当たるようになる(ただし上のずれは別途残る)。

★**2026-08-24 の未着手**: `pf_isOfIsotropicType`(`𝒞^pf` の根 1 の部分が
isotropic 型)はこの穴で止まっている。段取りは
`Found/FrdI/Prop53PfCatRoot.lean` の ★2 に書いてある。

### 14 の続報(2026-08-24、2 回目の試作)

**主因は「自前の同義語ほどき関数を書いたこと」**だった。
在庫に `pfDown (A : PfCat P F) : C := A` があるのに `pfObjDown` を自分で足すと、
`pfDiv` などの暗黙引数(在庫の `pfDown` で書かれている)と噛み合わない。

★**在庫の同義語ほどき関数に揃える**。揃えた上で

* `pfCat_comp_eq : f ≫ g = compPf P F f g := rfl`
* `toPfCat_map_eq : (toPfCat P F).map ψ = toHomPf ψ := rfl`

の 2 本を橋にすると、`pfDiv_comp` / `pfDeg_comp` / `rootDiv` / `rootBase` /
`compRoot_root_one` はすべて当たるようになった。

★★**目標の側で `rw` しない** —— `≫` を開いた目標は
`instances` 透明度で型が付かなくなる(エラーに
「The target expression is not type-correct under the `instances` transparency level」
と出る)。**自分で述べた `have` の側で `rw` して、最後に `exact` で defeq に頼る。**

★残るのは `IsIso` のインスタンス合成(既存の仮定が暗黙引数の書かれ方の違いで
拾われない)で、これは `haveI` を**目標に現れる形そのまま**で置き直すしかない。

### 14 の続報(2026-08-25、3 回目)—— `𝒟` 側の合成が `Functor.comp` をまたぐとき

`Cor54SeamCls.lean` で、**底 `𝒟` の射の合成**でも同じ穴に落ちた。
`sq : P₁.proj ⋙ ΨB ≅ Ψ ⋙ P₂.proj` の成分 `sq.hom.app X` の終域は
`(Ψ ⋙ P₂.proj).obj X` だが、`P₂.Base f` の始域は `(P₂.toElem.obj (Ψ.obj X)).base`。
**defeq だが `instances` 透明度では型が付かない**ので、

* `Category.assoc` / `IsIso.hom_inv_id` が `rw` で当たらない
* `IsIso (P₂.Base (Ψ.map φ))` が `haveI` で置いてあっても**合成に失敗する**

**対処(測定済み、3 点)**:

1. **`inv` は `@inv _ _ _ _ f h` でインスタンスを明示**する。
   `IsIso` は `Prop` なので、どの証明を渡しても defeq。合成に頼らないのが確実。
2. **`Category.assoc` / `hom_inv_id` は `have` で項として置く**
   (`(Category.assoc _ _ _).symm` / `@IsIso.hom_inv_id _ _ _ _ f h`)。
   置いた `have` を `rw` するか、`Eq.trans` / `congrArg` で項の側から組む。
3. **`rw` の末尾で閉じない**ことがある —— 目標が見た目 `X = X` になっても
   `rw` の自動 `rfl` は reducible 透明度なので通らない。**`rfl` を 1 行足す**。

★★**シェルの罠(同日)**: `perl -0pi -e 's/.../(@inv ...)/'` の置換文字列で
`@inv` が**配列展開されて消える**。置換に `@` を含めるときは `sed` を使うか
`\@` でエスケープする。Lean ファイルの中身は Write/Edit で書くのが安全。

## 14 の続報(2026-08-25、4 回目)—— `inv` の `rw` は**抽象補題へ逃がす**のが唯一安定

`Prop44Univ.lean`(`𝒞^birat` の普遍性)で 4 つ新しい失敗形に当たった。
**3 回目までの対処(`@inv` で明示、`have` で項にする)では足りない場面がある。**

| 失敗形 | 出るメッセージ | 直し方 |
|---|---|---|
| `rw [Category.assoc]` が当たらない | `Did not find an occurrence of the pattern (?f ≫ ?g) ≫ ?h` (**目に見えて在るのに**) | ↓ |
| `rw [IsIso.hom_inv_id]` が当たらない | `Did not find an occurrence of the pattern inv ?f ≫ ?f` | ↓ |
| `rw [Ω.map_comp]` / `rw [Ω.map_id]` | `motive is not type correct`(`IsIso _a` に型が付かない) | ↓ |

★**原因**は共通で、`IsIso` の引数がメタ変数だと**インスタンス探索が走らない**こと。
`rw` は補題 `IsIso.hom_inv_id` を使うのに `[IsIso ?f]` を解決できず、
パターンそのものを作れないので「見つからない」と言う。

### ★★対処 —— **圏の中だけの抽象補題**に逃がす

`IsIso` を**インスタンス束縛**にした補題を別に立てて、`exact` で当てる。
補題の中では `?f` が本物の変数なので `rw` が普通に動く。

```lean
theorem frac_key_aux {X Y T U : E} (g : X ⟶ Y) [IsIso g] (a : X ⟶ T) (p : Y ⟶ U)
    (w : T ⟶ U) [IsIso w] (hsq : g ≫ p = a ≫ w) : inv g ≫ a = p ≫ inv w := by
  rw [IsIso.inv_comp_eq, ← Category.assoc, hsq, Category.assoc, IsIso.hom_inv_id,
    Category.comp_id]
```

**`rw [Ω.map_comp]` の motive 問題**も、等式を仮定に外出しすれば消える:

```lean
theorem frac_comp_aux ... (gz : X ⟶ Z) [IsIso gz] (aq : X ⟶ V)
    (hgz : gz = g ≫ z) (haq : aq = a ≫ q) ... := by
  subst hgz; subst haq   -- ★ここで初めて合成の形になる。inv の下を rw しない
```

### ★注意 1 —— 抽象補題の **universe を 1 本に固定する**

`D` / `C` を使わない補題を同じ section に置くと、`Category.{max v u2 v2} E` の
`max v u2 v2` が**そのまま universe 変数 3 本として一般化**され、
呼び出し側で `stuck at solving universe constraint` になる。

```lean
section FracAux
universe uv
variable {E : Type uE} [Category.{uv} E]   -- ★ 1 本にする
```

### ★注意 2 —— 呼ぶときは **インスタンス引数を `@` で明示**する

`exact frac_comp_aux g z a p w q …` は
`failed to synthesize IsIso (Ω.map …)` で落ちる。
**局所インスタンスは在る**(`have test : IsIso … := inferInstance` は通る)のに、
暗黙の**対象**がメタ変数のうちに探索が走るため。

```lean
exact @frac_comp_aux _ _ _ _ _ _ _ _ g hγ z hZ a p w hW q gz hP aq
  (Ω.map_comp _ _) (Ω.map_comp _ _) key
```

## 15. `calc` は defeq を渡らない —— `have` ＋ `Eq.trans` に置き換える

`Prop55BiratOmega.lean` で繰り返し当たった形。

```
error: invalid 'calc' step, failed to synthesize `Trans` instance
  Trans Eq Eq ?m.848
```

★**原因**: `calc` の各段は**型が構文的に一致**していないと `Trans` が付かない。
`omegaObj F F' ⟨A,k⟩` と `⟨biratUp A, k⟩` のように **defeq だが構文が違う**と落ちる。

★**逃げ方**: 段を `have s1 := …` で置き(型注釈を書かない)、
最後に `(((s1.trans s2).trans s3) … )` で繋ぐ。
`Eq.trans` は defeq で通るので、これだけで直る。

```lean
  have s1 := congrArg (compRoot Q F' _) (omegaMap_pfKappa F F' B k).symm
  have s2 := (omegaMap_comp F F' _ _).symm
  have s3 := congrArg (omegaMap F F') (rootMap_spec (F := F) hfi f k)
  exact rootMap_ext (F := F') hfiB _ _ (((s1.trans s2).trans s3).trans hR.symm)
```

## 16. `IsIso.inv_eq_of_hom_inv_id` は `@` で `f` とインスタンスを明示する

```
error: failed to synthesize instance of type class IsIso (rtExt (biratPre P G) F' (biratUp P G A) 1)
```

局所インスタンスは在るのに落ちる。★**暗黙の `f` がメタ変数のうちに探索が走る**ため
(idiom 14 の注意 2 と同じ)。

★**引数の順は `{f} [IsIso f] {g}`** —— インスタンスは `g` より**前**である。

```lean
refine (@IsIso.inv_eq_of_hom_inv_id _ _ _ _
  (rtExt (biratPre P G) F' (biratUp P G A) 1) hq1   -- ★f, インスタンスの順
  (biratRtIso F F' A 1 ≫ (toBiratCat P G).map (rtOneInv A)) ?_).symm
```

## 17. `set_option … in` は docstring の**前**に置く

`variable (F) in` と併用するときは

```lean
set_option maxHeartbeats 1600000 in
variable (F) in
/-- doc -/
theorem foo …
```

★docstring の中に紛れ込ませると本文の一部になって効かない
(`sed` で行番号を数えて挿入すると起きやすい)。

## 18. `rw` は「対象が関手の像の形」だと当たらない —— 項スタイルへ

`Prop55BiratOmega.lean` の組み立てで**丸一日ぶんの試行を溶かした形**。測定結果:

補題の側の対象が `omegaObj F F' X`(関手の対象写像を当てた形)で、
目標の側では既に `⟨biratUp Z₀, k⟩` に簡約されているとき、`rw` は

```
omegaObj F F' ?X =?= ⟨biratUp Z₀, k⟩
```

を解こうとする。`omegaObj` は構造体リテラルを返す `def` なので、
これは `?X.obj =?= Z₀` という**メタ変数の射影**に化けて解けない。
★エラーは「`Did not find an occurrence of the pattern`」と出るが、**原因はこれ**。

★★**対処**: 組み立ては `Eq.trans` / `congrArg` の**項スタイル**で書く。
`≫` 版の補題を用意しても同じ理由で `rw` は当たらない(実測)。

### ★付随して測ったこと

* `HomRoot P F X Y` は `≫` の**左**に置けば `Quiver.Hom` に解けるが、
  等式の**右辺**に単独で置くと解けない(メタ変数の postpone)。
  `≫` 版の補題を書くときは型注釈
  `(… : omegaObj F F' X ⟶ omegaObj F F' Zz)` が要る。
* `Mono ((Ω).mapIso i).hom` の**インスタンス探索は heartbeat を使い切る**。
  `Iso.cancel_iso_hom_right` を使うこと。
* 文の中で `show T from e` を使うと `have this := e; this` に展開され、
  以後の `Category.assoc` などの単一化を**全部止める**。
  文では名前つき引数 `(X := …) (Y := …)` を使い、`show` は**戦術の側**に置く。

## `letI := algOfHom f` で `Algebra L L` を入れても `Algebra (𝓞 L) (𝓞 L)` は f を見ない

**失敗形**

```lean
letI := algOfHom f   -- Algebra L.toIF L.toIF
have h : ((algebraMap (𝓞 L.toIF) (𝓞 L.toIF)) x : L.toIF) = (FinSub.hom f) (x : L.toIF) := rfl
-- Type mismatch: rfl has type ?m = ?m
```

★★原因: `Algebra (𝓞 A) (𝓞 A)` には **`Algebra.id` が既にインスタンスとしてある**ので、
`letI` で入れた `Algebra L L` から**派生する**インスタンスより先に選ばれる。
その結果 `algebraMap (𝓞 L) (𝓞 L)` は**恒等写像**に解決され、`f` を一切見ない。

★★★**ただし既存の定義（`resHOS` など）は無事**である ——
`resHOS {L M : Type} [Algebra L M]` の本体は `L`・`M` が**相異なる型変数**の状態で
エラボレートされており、そこでは `Algebra.id` は候補にならないので
**派生インスタンスが焼き付いている**。あとから `L := M := L.toIF` を代入しても
型クラス探索は**やり直されない**ので、`resHOS` は正しく `f` に依存する。

**対処**: `L = M` の場合を扱う補題は、**一般の `{L M} [Algebra L M]` の側に置く**。

```lean
theorem asIdeal_resHOS (V) :
    (resHOS (L := L) V).asIdeal = V.asIdeal.comap (algebraMap (𝓞 L) (𝓞 M)) := rfl
```

そこでは `((algebraMap (𝓞 L) (𝓞 M)) x : M) = algebraMap L M (x : L)` が `rfl` である
（mathlib の `RingOfIntegers.instAlgebra` はそう作られている）。
★同じ式を `L = M` に代入した文脈で**書き直す**と `Algebra.id` を拾って壊れる。

## `bernoulli` の値は `decide` でも `norm_num [bernoulli]` でも出ない

**失敗形**:

```lean
example : (bernoulli 4 : ℚ) = -1/30 := by decide
-- Decidable インスタンスが `bernoulli' 4` の再帰で止まる
example : (bernoulli 6 : ℚ) = 1/42 := by norm_num [bernoulli]
-- `bernoulli.eq_1` と `bernoulli'_def` が looping simp theorem になり maxRecDepth
```

★原因: `bernoulli'` は `bernoulli' n = 1 - ∑_{k<n} C(n,k)/(n-k+1) * bernoulli' k` という
**自分自身を含む** well-founded 再帰。simp に渡すと展開が止まらない。

**対処**: `bernoulli'_def` を **1 回だけ `rw`** して、下の値は**既に証明した補題**として
`norm_num` に渡す。`Nat.choose` も明示的に渡さないと `Nat.choose 4 2` が残る。

```lean
theorem bernoulli'_four : (bernoulli' 4 : ℚ) = -1/30 := by
  rw [bernoulli'_def]
  norm_num [Finset.sum_range_succ, bernoulli'_zero, bernoulli'_one, bernoulli'_two,
    bernoulli'_three, Nat.choose]
```

`n = 3, 4, 5, 6` と 1 段ずつ積み上げる。`bernoulli n` へは
`bernoulli_eq_bernoulli'_of_ne_one`(`n ≠ 1`)で移る。
★`riemannZeta 6` も mathlib に無いので、`riemannZeta_two_mul_nat` と `B₆` から作る
(`riemannZeta_two`・`riemannZeta_four` はある)。

## `Ring.inverse` は環準同型と交換しない(単元でない限り)

**失敗形**: `Ideal.Quotient.mk` や特殊化 `φ : R →+* R'` を `tateXterm t = t * Ring.inverse (1-t)^2`
のような式に当てて、`φ (Ring.inverse x) = Ring.inverse (φ x)` を暗黙に使ってしまう。

★`Ring.inverse x` は `x` が単元でないとき**既定値 `0`** を返す。`φ x` が単元になっても
`φ 0 = 0 ≠ Ring.inverse (φ x)` なので**一般には成り立たない**。

**対処**: `IsUnit x` を仮定に持つ。

```lean
theorem map_ring_inverse (φ : R →+* R') {x : R} (hx : IsUnit x) :
    φ (Ring.inverse x) = Ring.inverse (φ x) := by
  have hux : IsUnit (φ x) := hx.map φ
  have h : φ (Ring.inverse x) * φ x = 1 := by
    rw [← map_mul, Ring.inverse_mul_cancel _ hx, map_one]
  calc φ (Ring.inverse x) = φ (Ring.inverse x) * (φ x * Ring.inverse (φ x)) := by
        rw [Ring.mul_inverse_cancel _ hux, mul_one]
    _ = (φ (Ring.inverse x) * φ x) * Ring.inverse (φ x) := by ring
    _ = Ring.inverse (φ x) := by rw [h, one_mul]
```

★`Ring.eq_inverse_of_mul_eq_one_left` は**無い**ので上の `calc` で作る。
★★これが理由で、`Ring.inverse` を含む式の mod `I` 議論では `Ideal.Quotient` を使わず
**差の分解**(`Y²−g² = (Y−g)(Y+g)` 等)で処理する。
